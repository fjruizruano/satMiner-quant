#!/usr/bin/env python3
"""
satminer_quant_py3_allinone.py

Python 3 port + "single-script" integration of the satMiner quantification protocol.

This script integrates (and modernizes) the following python2 scripts from:
  - https://github.com/fjruizruano/ngs-protocols

Integrated functionality:
  - divsum_ab.py        (now implemented as divsum_ab())
  - sat_subfam2fam.py   (now implemented as sat_subfam2fam(); still calls calcDivergenceFromAlign.pl)
  - divsum_to_rl.py     (now implemented as divsum_to_rl(); still calls Rscript for plotting)
  - replace_patterns.py (now implemented as replace_patterns())

Notes / external dependencies that remain:
  - calcDivergenceFromAlign.pl (Perl, from RepeatMasker utilities) is still invoked.
  - R + ggplot2, reshape2, plyr, RColorBrewer are still invoked (same as original pipeline).
  - Biopython is required (Bio.SeqIO), same as original.

Usage (compatible with the original):
    satminer_quant_py3_allinone.py SamplesFile FastaFileMonomers

The SamplesFile format is the same as the original script:
  First line:  sp_name<TAB>ref_library<TAB>rep_land
  Next lines:  library_name<TAB>divsum_file<TAB>nucs

rep_land can be "NO" or a comma-separated list of pairs like: lib1-lib2,lib3-lib4
"""

from __future__ import annotations

import argparse
import operator
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

try:
    from Bio import SeqIO
except ImportError as e:
    raise SystemExit(
        "Biopython is required (pip install biopython). "
        "The original pipeline uses Bio.SeqIO as well."
    ) from e


# -------------------------
# helpers
# -------------------------

def run(cmd: List[str] | str, *, shell: bool = False) -> None:
    """Run a command (prints it first) and raise if it fails."""
    if isinstance(cmd, list):
        printable = " ".join(cmd)
    else:
        printable = cmd
    print(printable)
    subprocess.run(cmd, shell=shell, check=True)


def read_lines(path: Path) -> List[str]:
    return path.read_text(encoding="utf-8", errors="replace").splitlines(True)


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def parse_patterns(pattern_file: Path) -> Dict[str, str]:
    patterns: Dict[str, str] = {}
    for line in read_lines(pattern_file):
        parts = line.split()
        if len(parts) >= 2:
            patterns[parts[0]] = parts[1]
    return patterns


def fasta_id_to_length(fasta_path: Path) -> Dict[str, int]:
    """Return a mapping {record_id: sequence_length} for the monomer FASTA."""
    lengths: Dict[str, int] = {}
    if not fasta_path.exists():
        return lengths
    # Use Biopython if available (it is required by this pipeline).
    with fasta_path.open("r", encoding="utf-8", errors="replace") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            # SeqIO record.id is up to first whitespace, which matches typical usage.
            lengths[str(rec.id)] = len(rec.seq)
    return lengths


def load_table_lengths(table_path: Path) -> Dict[str, str]:
    """Read table.txt-like files as {record_id: monomer_length_as_string}."""
    lengths: Dict[str, str] = {}
    if not table_path.exists():
        return lengths
    with table_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                lengths[parts[0]] = parts[1]
    return lengths


def split_alias_core(alias_id: str) -> str:
    """
    Remove the trailing '-<digits>' length suffix from an intermediate alias id.
    Example: 'DvittaR2CL250b_A1B-105' -> 'DvittaR2CL250b_A1B'
    """
    if "-" not in alias_id:
        return alias_id
    left, right = alias_id.rsplit("-", 1)
    return left if right.isdigit() else alias_id


def match_family_from_alias(alias_core: str, family_ids_sorted: List[str]) -> Tuple[str | None, str]:
    """
    Match an intermediate alias such as 'LeaderA' or 'LeaderB' back to its family leader id.
    Returns (leader_id, variant_suffix_letters).
    """
    for fam_id in family_ids_sorted:
        if alias_core == fam_id:
            return fam_id, ""
        if alias_core.startswith(fam_id):
            suffix = alias_core[len(fam_id):]
            if suffix.isalpha():
                return fam_id, suffix
    return None, ""




def cluster_key_from_id(seq_id: str) -> str | None:
    """Return a cluster key for multi-variant monomer ids like F1V2..._A1.

    Examples:
      F1V2dvR5CL397b_A1 -> F1_A1
      F14V4dvFR2CL249_B1 -> F14_B1
    Returns None for ids that do not match this multi-variant pattern.
    """
    m = re.match(r'^(F\d+)V\d+.*_([A-Z]\d+)$', seq_id)
    if not m:
        return None
    return f"{m.group(1)}_{m.group(2)}"


def rebuild_pattern_from_original_abundances(abund_path: Path, fasta_monomers: Path) -> None:
    """Rebuild pattern.txt/selection.txt/table.txt using all original FASTA variants.

    Strategy:
      - group multi-variant ids using the original FASTA id pattern (e.g. F1V1..._A1, F1V2..._A1)
      - rank members by decreasing abundance in the primary library; if all are 0, use the secondary
      - write pattern.txt with every non-leader member -> leader, even if the member is absent from the primary
    This ensures secondary-only variants are still included later when editing *.align files.
    """
    if not abund_path.exists():
        print(f"WARNING: {abund_path} not found; keeping the existing pattern.txt")
        return

    rows = []
    with abund_path.open("r", encoding="utf-8", errors="replace") as fh:
        _header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rid = parts[0]
            try:
                length = int(parts[1])
            except Exception:
                length = 0
            try:
                primary_ab = float(parts[2])
            except Exception:
                primary_ab = 0.0
            try:
                secondary_ab = float(parts[4])
            except Exception:
                secondary_ab = 0.0
            rows.append((rid, length, primary_ab, secondary_ab))

    fasta_records = list(SeqIO.parse(str(fasta_monomers), "fasta"))
    fasta_by_id = {str(rec.id): rec for rec in fasta_records}
    ids_in_fasta = [str(rec.id) for rec in fasta_records]

    row_by_id = {rid: (length, pab, sab) for rid, length, pab, sab in rows}
    groups: Dict[str, List[str]] = {}
    for rid in ids_in_fasta:
        key = cluster_key_from_id(rid) or rid
        groups.setdefault(key, []).append(rid)

    def rank_key(rid: str) -> tuple[int, float, float, str]:
        _length, pab, sab = row_by_id.get(rid, (len(str(fasta_by_id[rid].seq)), 0.0, 0.0))
        if pab > 0:
            tier = 0
        elif sab > 0:
            tier = 1
        else:
            tier = 2
        return (tier, -pab, -sab, rid)

    pattern_lines: List[str] = []
    selection_ids: List[str] = []
    var_number: Dict[str, int] = {}

    for _key, members in groups.items():
        if len(members) == 1 and cluster_key_from_id(members[0]) is None:
            rid = members[0]
            selection_ids.append(rid)
            var_number[rid] = 1
            continue

        ordered = sorted(members, key=rank_key)
        leader = ordered[0]
        selection_ids.append(leader)
        var_number[leader] = len(members)
        for member in ordered[1:]:
            pattern_lines.append(f"{member}	{leader}")

    Path("pattern.txt").write_text("\n".join(pattern_lines) + ("\n" if pattern_lines else ""), encoding="utf-8")
    Path("selection.txt").write_text("\n".join(selection_ids) + ("\n" if selection_ids else ""), encoding="utf-8")

    selection_extract_path = Path("selection.txt.extract")
    with selection_extract_path.open("w", encoding="utf-8") as out:
        for rid in selection_ids:
            rec = fasta_by_id.get(rid)
            if rec is not None:
                SeqIO.write(rec, out, "fasta")

    with Path("table.txt").open("w", encoding="utf-8") as out:
        for rid in selection_ids:
            rec = fasta_by_id[rid]
            seq = str(rec.seq)
            length = len(seq)
            at = seq.count("A") + seq.count("T") + seq.count("a") + seq.count("t")
            at_perc = float(at / length) if length else 0.0
            out.write(f"{rid}\t{length}\t{at_perc}\t{var_number.get(rid, 1)}\n")

def parse_divsum_file(divsum_path: Path, nucs: int) -> dict[str, object]:
    """Parse one RepeatMasker .divsum file into divergence, absolute counts, and relative abundances."""
    data = read_lines(divsum_path)
    s_matrix_div = data.index("-----\t------\t------\t-----------\t-------\n")
    s_matrix = data.index("Coverage for each repeat class and divergence (Kimura)\n")

    divergence: Dict[str, str] = {}
    for line in data[s_matrix_div + 1 : s_matrix - 2]:
        info = line.split()
        if len(info) >= 2:
            divergence[info[1]] = info[-1]

    elements = data[s_matrix + 1].split()
    matrix_mut = [[element.split("/")[-1], []] for element in elements[1:]]
    n_el = len(matrix_mut)
    for line in data[s_matrix + 2 : s_matrix + 45]:
        info = line.split()[1:]
        if len(info) < n_el:
            continue
        for n in range(n_el):
            matrix_mut[n][1].append(int(info[n]))

    family_abs = {matrix_mut[n][0]: sum(matrix_mut[n][1]) for n in range(n_el)}
    family_rel = {k: (v / nucs) for k, v in family_abs.items()}
    return {
        "nucs": nucs,
        "divergence": divergence,
        "family_abs": family_abs,
        "family_rel": family_rel,
    }


def parse_original_divsum_counts(samples_file: Path) -> tuple[str, list[str], dict[str, dict[str, object]]]:
    """Parse the ORIGINAL .divsum files (not .align.fam.divsum) for all libraries."""
    samples = read_lines(samples_file)
    samples_rl = samples[0].rstrip("\n").split("\t")
    ref_library = samples_rl[1]

    lib_dict: Dict[str, List[str]] = {}
    for lib in samples[1:]:
        lib = lib.rstrip("\n").split("\t")
        if not lib or len(lib) < 3:
            continue
        lib_dict[lib[0]] = lib[1:]

    parsed: dict[str, dict[str, object]] = {}
    lib_order = list(lib_dict.keys())
    for library, vals in lib_dict.items():
        divsum_path = Path(vals[0])
        nucs = int(vals[1])
        parsed[library] = parse_divsum_file(divsum_path, nucs)

    return ref_library, lib_order, parsed

def load_pattern_leaders(pattern_file: Path, fasta_ids: Iterable[str]) -> Dict[str, str]:
    """Return {original_id: cluster_leader_id} using pattern.txt member->leader mappings."""
    leaders: Dict[str, str] = {rid: rid for rid in fasta_ids}
    if not pattern_file.exists():
        return leaders
    direct: Dict[str, str] = {}
    for line in read_lines(pattern_file):
        parts = line.split()
        if len(parts) >= 2:
            direct[parts[0]] = parts[1]

    def resolve(rid: str) -> str:
        seen = set()
        cur = rid
        while cur in direct and cur not in seen:
            seen.add(cur)
            cur = direct[cur]
        return cur

    for rid in list(leaders):
        leaders[rid] = resolve(rid)
    return leaders


def write_final_renamed_fasta(
    source_fasta: Path,
    out_fasta: Path,
    eq_map: Dict[str, str],
    leader_len_map: Dict[str, str],
    source_len_map: Dict[str, int],
    leader_by_src: Dict[str, str],
    ref_counts: Dict[str, int],
    secondary_counts: Dict[str, int],
    ref_rel: Dict[str, float],
    secondary_rel: Dict[str, float],
    ref_divergence: Dict[str, str],
    secondary_divergence: Dict[str, str],
    equivalence_out: Path | None = None,
    abundance_out: Path | None = None,
    abundance_abs_out: Path | None = None,
    ref_library: str | None = None,
    secondary_library: str | None = None,
) -> None:
    """
    Build final renamed FASTA using all original FASTA variants plus abundances from the ORIGINAL divsums.
    """
    if not source_fasta.exists():
        print(f"WARNING: source FASTA not found: {source_fasta}")
        return

    source_records = list(SeqIO.parse(str(source_fasta), "fasta"))
    source_ids = [str(rec.id) for rec in source_records]

    base_by_src: Dict[str, str] = {}
    for src_id in source_ids:
        leader_id = leader_by_src.get(src_id, src_id)
        if src_id in eq_map:
            base_by_src[src_id] = eq_map[src_id]
        elif leader_id in eq_map:
            base_by_src[src_id] = eq_map[leader_id]
        else:
            base_by_src[src_id] = src_id

    members_by_leader: Dict[str, List[str]] = {}
    for src_id in source_ids:
        members_by_leader.setdefault(leader_by_src.get(src_id, src_id), []).append(src_id)

    suffix_by_src: Dict[str, str] = {}
    for leader_id, members in members_by_leader.items():
        leader_base = eq_map.get(leader_id, base_by_src.get(leader_id, leader_id))
        keep_under_leader = [m for m in members if base_by_src[m] == leader_base]

        def rank_key(src_id: str) -> tuple[int, float, float, str]:
            ref_ab = ref_counts.get(src_id, 0)
            sec_ab = secondary_counts.get(src_id, 0)
            if ref_ab > 0:
                tier = 0
            elif sec_ab > 0:
                tier = 1
            else:
                tier = 2
            return (tier, -ref_ab, -sec_ab, src_id)

        ordered_keep = sorted(keep_under_leader, key=rank_key)
        if len(members) == 1:
            suffix_by_src[members[0]] = ""
        elif len(ordered_keep) == 1:
            suffix_by_src[ordered_keep[0]] = "A"
        else:
            for idx, src_id in enumerate(ordered_keep):
                suffix_by_src[src_id] = suffix_from_index(idx)

        for m in members:
            if m not in suffix_by_src:
                suffix_by_src[m] = ""

    equivalence_pairs: List[Tuple[str, str]] = []
    renamed_records: List[Tuple[str, str, str]] = []
    abundance_rows: List[Tuple[str, str, float, str, float, str]] = []
    abundance_abs_rows: List[Tuple[str, str, int, str, int, str]] = []

    for src_rec in source_records:
        src_id = str(src_rec.id)
        seq = str(src_rec.seq)
        variant_len = len(seq)
        final_base = base_by_src[src_id]
        suffix_for_name = suffix_by_src.get(src_id, "")
        leader_id = leader_by_src.get(src_id, src_id)

        leader_len_id = src_id if src_id in eq_map else leader_id
        if leader_len_id in leader_len_map:
            leader_len = leader_len_map[leader_len_id]
        elif leader_len_id in source_len_map:
            leader_len = str(source_len_map[leader_len_id])
        else:
            leader_len = str(variant_len)

        variant_name = f"{final_base}{suffix_for_name}-{variant_len}"
        family_name = f"{final_base}-{leader_len}"
        equivalence_pairs.append((src_id, variant_name))
        renamed_records.append((variant_name, family_name, seq))
        abundance_rows.append((src_id, variant_name, float(ref_rel.get(src_id, 0.0)), str(ref_divergence.get(src_id, 'NA')), float(secondary_rel.get(src_id, 0.0)), str(secondary_divergence.get(src_id, 'NA'))))
        abundance_abs_rows.append((src_id, variant_name, int(ref_counts.get(src_id, 0)), str(ref_divergence.get(src_id, 'NA')), int(secondary_counts.get(src_id, 0)), str(secondary_divergence.get(src_id, 'NA'))))

    renamed_records.sort(key=lambda x: x[0])
    equivalence_pairs.sort(key=lambda x: x[1])
    abundance_rows.sort(key=lambda x: x[1])
    abundance_abs_rows.sort(key=lambda x: x[1])

    with out_fasta.open("w", encoding="utf-8") as out_handle:
        for variant_name, family_name, seq in renamed_records:
            out_handle.write(f">{variant_name}#Satellite/{family_name}\n{seq}\n")

    if equivalence_out is not None:
        with equivalence_out.open("w", encoding="utf-8") as out_eq:
            for old_id, new_id in equivalence_pairs:
                out_eq.write(f"{old_id}\t{new_id}\n")

    ref_name = ref_library or "primary"
    sec_name = secondary_library or "secondary"
    if abundance_out is not None:
        with abundance_out.open("w", encoding="utf-8") as out_ab:
            out_ab.write(f"OriginalID\tFinalID\t{ref_name}_abundance\t{ref_name}_divergence\t{sec_name}_abundance\t{sec_name}_divergence\n")
            for old_id, final_id, ref_ab, ref_div, sec_ab, sec_div in abundance_rows:
                out_ab.write(f"{old_id}\t{final_id}\t{ref_ab}\t{ref_div}\t{sec_ab}\t{sec_div}\n")

    if abundance_abs_out is not None:
        with abundance_abs_out.open("w", encoding="utf-8") as out_ab_abs:
            out_ab_abs.write(f"OriginalID\tFinalID\t{ref_name}_abundance\t{ref_name}_divergence\t{sec_name}_abundance\t{sec_name}_divergence\n")
            for old_id, final_id, ref_ab, ref_div, sec_ab, sec_div in abundance_abs_rows:
                out_ab_abs.write(f"{old_id}\t{final_id}\t{ref_ab}\t{ref_div}\t{sec_ab}\t{sec_div}\n")

def rename_fasta_and_dim_outputs(fasta_monomers: Path, samples_file: Path) -> None:
    """Create the final renamed FASTA outputs and per-sequence abundance tables using the ORIGINAL divsums."""
    eq_map = parse_patterns(Path("equivalences.txt"))
    if not eq_map:
        print("WARNING: equivalences.txt is empty or missing; final FASTA renaming skipped.")
        return

    leader_len_map = load_table_lengths(Path("table.txt"))
    source_len_map = fasta_id_to_length(fasta_monomers)
    ref_library, lib_order, parsed_original = parse_original_divsum_counts(samples_file)

    secondary_library = next((lib for lib in lib_order if lib != ref_library), None)
    ref_counts = parsed_original[ref_library]["family_abs"] if ref_library in parsed_original else {}
    secondary_counts = parsed_original[secondary_library]["family_abs"] if secondary_library and secondary_library in parsed_original else {}
    ref_rel = parsed_original[ref_library]["family_rel"] if ref_library in parsed_original else {}
    secondary_rel = parsed_original[secondary_library]["family_rel"] if secondary_library and secondary_library in parsed_original else {}
    ref_divergence = parsed_original[ref_library]["divergence"] if ref_library in parsed_original else {}
    secondary_divergence = parsed_original[secondary_library]["divergence"] if secondary_library and secondary_library in parsed_original else {}

    fasta_ids = [str(rec.id) for rec in SeqIO.parse(str(fasta_monomers), "fasta")]
    leader_by_src = load_pattern_leaders(Path("pattern.txt"), fasta_ids)

    final_fasta = fasta_monomers.with_name(fasta_monomers.name + ".fam")
    write_final_renamed_fasta(fasta_monomers, final_fasta, eq_map, leader_len_map, source_len_map, leader_by_src, ref_counts, secondary_counts, ref_rel, secondary_rel, ref_divergence, secondary_divergence, Path("equivalences.txt.fam"), Path("final_variant_abundances.txt"), Path("final_variant_abundances_absolute.txt"), ref_library, secondary_library)

    dim_fasta = fasta_monomers.with_name(fasta_monomers.name + ".dim")
    final_dim_fasta = fasta_monomers.with_name(fasta_monomers.name + ".dim.fam")
    if dim_fasta.exists():
        write_final_renamed_fasta(dim_fasta, final_dim_fasta, eq_map, leader_len_map, source_len_map, leader_by_src, ref_counts, secondary_counts, ref_rel, secondary_rel, ref_divergence, secondary_divergence, None, None, None, ref_library, secondary_library)
    else:
        print(f"WARNING: expected {dim_fasta} to exist for final .dim.fam output, but it was not found.")



def suffix_from_index(idx: int) -> str:
    """Return Excel-like variant suffixes: 0->A, 1->B, ..., 25->Z, 26->AA, ..."""
    if idx < 0:
        raise ValueError("idx must be >= 0")
    out = []
    n = idx
    while True:
        n, rem = divmod(n, 26)
        out.append(chr(ord('A') + rem))
        if n == 0:
            break
        n -= 1
    return ''.join(reversed(out))

def replace_patterns(input_file: Path, pattern_file: Path, *, output_suffix: str = ".fam") -> Path:
    """Replicates replace_patterns.py behavior: naive .replace for each key across each line."""
    patterns = parse_patterns(pattern_file)
    out_path = input_file.with_name(input_file.name + output_suffix)
    with input_file.open("r", encoding="utf-8", errors="replace") as r, out_path.open(
        "w", encoding="utf-8"
    ) as w:
        for line in r:
            for old, new in patterns.items():
                line = line.replace(old, new)
            w.write(line)
    return out_path


def extract_sequences(fasta_path: Path, ids: Iterable[str], out_path: Path) -> None:
    """Replacement for extract_seq.py: write FASTA records whose id is in 'ids'."""
    wanted = set(ids)
    with out_path.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            if rec.id in wanted:
                SeqIO.write(rec, out, "fasta")


# -------------------------
# integrated scripts
# -------------------------

def divsum_ab(divsum_file: Path, fasta_monomers: Path) -> None:
    """
    Port of ngs-protocols/divsum_ab.py.

    Produces (same filenames as original, in current working directory):
      - {divsum_file}.counts
      - pattern.txt
      - selection.txt
      - selection.txt.extract
      - table.txt
      - {fasta_monomers}.abc
      - {fasta_monomers}.dim.abc   (expects {fasta_monomers}.dim to exist)

    It also implicitly defines "families" by grouping repeat names whose last underscore
    token has length 2, selecting the most abundant as the leader.
    """
    data = read_lines(divsum_file)

    header = "Coverage for each repeat class and divergence (Kimura)\n"
    try:
        matrix_start = data.index(header)
    except ValueError as e:
        raise ValueError(
            f"Could not find expected header line in {divsum_file}: {header.strip()}"
        ) from e

    matrix = data[matrix_start + 1 :]

    names_line = matrix[0]
    info = names_line.split()
    fams = info[1:]  # skip first column

    li_clean: List[List[str]] = []
    for fam in fams:
        fam_clean = fam.split("/")
        # original script uses fam_clean[1]
        li_clean.append([fam_clean[1] if len(fam_clean) > 1 else fam_clean[0]])

    info_len = len(li_clean)

    for line in matrix[1:]:
        parts = line.split()
        row = parts[1:]
        if len(row) < info_len:
            continue
        for i in range(info_len):
            li_clean[i].append(row[i])

    # write .counts and build counts list
    counts: List[Tuple[str, int]] = []
    counts_path = divsum_file.with_name(divsum_file.name + ".counts")
    with counts_path.open("w", encoding="utf-8") as out:
        for el in li_clean:
            name = el[0]
            numbers = [int(x) for x in el[1:] if x.strip() != ""]
            total = sum(numbers)
            counts.append((name, total))
            out.write(f"{name}\t{total}\n")

    fam_count: Dict[str, Tuple[str, int]] = {}
    fam_all: Dict[str, Dict[str, int]] = {}
    fam_members: Dict[str, List[str]] = {}
    seq_list: List[str] = []
    var_number: Dict[str, int] = {}

    for clase, count in counts:
        clase_fam = clase.split("_")[-1]
        if len(clase_fam) == 2:
            if clase_fam in fam_count:
                fam_members[clase_fam].append(clase)
                last_leader, last_count = fam_count[clase_fam]
                if count > last_count:
                    fam_count[clase_fam] = (clase, count)
                fam_all[clase_fam][clase] = count
            else:
                fam_members[clase_fam] = [clase]
                fam_count[clase_fam] = (clase, count)
                fam_all[clase_fam] = {clase: count}
        else:
            seq_list.append(clase)
            var_number[clase] = 1

    # pattern.txt
    pattern_path = Path("pattern.txt")
    with pattern_path.open("w", encoding="utf-8") as out:
        for fam_key in fam_count:
            leader = fam_count[fam_key][0]
            seq_list.append(leader)
            members = fam_members[fam_key]
            var_number[leader] = len(members)
            for member in members:
                if member != leader:
                    out.write(f"{member}\t{leader}\n")

    # selection.txt
    selection_path = Path("selection.txt")
    write_text(selection_path, "\n".join(seq_list))

    # selection.txt.extract (replacement for extract_seq.py)
    selection_extract_path = Path("selection.txt.extract")
    extract_sequences(fasta_monomers, seq_list, selection_extract_path)

    # table.txt
    table_path = Path("table.txt")
    with table_path.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(selection_extract_path), "fasta"):
            seq = str(rec.seq)
            rid = str(rec.id)
            length = len(seq)
            at = seq.count("A") + seq.count("T") + seq.count("a") + seq.count("t")
            at_perc = float(at / length) if length else 0.0
            v_number = var_number.get(rid, 1)
            out.write(f"{rid}\t{length}\t{at_perc}\t{v_number}\n")

    # abc naming
    abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    abc_dict: Dict[str, str] = {}
    for fam_key, family in fam_all.items():
        # python3: items()
        family_sort = sorted(family.items(), key=operator.itemgetter(1), reverse=True)
        for i, (var, _) in enumerate(family_sort):
            if i >= len(abc):
                # keep behavior simple: stop if more than 26 variants
                break
            leader = fam_count[fam_key][0]
            abc_dict[var] = leader + abc[i]

    # write fasta.abc
    abc_path = fasta_monomers.with_name(fasta_monomers.name + ".abc")
    len_dict: Dict[str, int] = {}
    with abc_path.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(fasta_monomers), "fasta"):
            new = rec.id
            length = len(str(rec.seq))
            if rec.id in abc_dict:
                new = abc_dict[rec.id]
            out.write(f">{new}-{length}\n{rec.seq}\n")
            len_dict[new] = length

    # write fasta.dim.abc (expects fasta.dim)
    dim_path = fasta_monomers.with_name(fasta_monomers.name + ".dim")
    dim_abc_path = fasta_monomers.with_name(fasta_monomers.name + ".dim.abc")
    if dim_path.exists():
        with dim_abc_path.open("w", encoding="utf-8") as out:
            for rec in SeqIO.parse(str(dim_path), "fasta"):
                new = rec.id
                if rec.id in abc_dict:
                    new = abc_dict[rec.id]
                length = len_dict.get(new, len(str(rec.seq)))
                out.write(f">{new}-{length}\n{rec.seq}\n")
    else:
        print(f"WARNING: expected {dim_path} to exist for .dim.abc output, but it was not found.")


def sat_subfam2fam(align_file: Path, pattern_file: Path) -> Path:
    """
    Port of ngs-protocols/sat_subfam2fam.py.

    Creates:
      - {align_file}.fam
      - {align_file}.fam.divsum (via calcDivergenceFromAlign.pl)

    Returns the created .fam file path.
    """
    fam_path = replace_patterns(align_file, pattern_file, output_suffix=".fam")
    # keep same external call as original
    run(f"calcDivergenceFromAlign.pl -s {align_file}.fam.divsum {fam_path}", shell=True)
    return fam_path


def divsum_to_rl(samples_file: Path, fasta_monomers: Path) -> None:
    """
    Port of ngs-protocols/divsum_to_rl.py.

    Reads table.txt (if present) and writes:
      - equivalences.txt
      - <library>.abdiv, <library>_rl.txt, <library>_rl.R, <library>_rl.pdf (if Rscript works)
      - <rep_land_pair>_rl.* for subtractive landscapes (if configured)

    This stays very close to the original behavior (including generating and running R scripts).
    """
    # loading samples file
    samples = read_lines(samples_file)
    samples_rl = samples[0].rstrip("\n").split("\t")
    sp_name = samples_rl[0]
    ref_library = samples_rl[1]
    rep_land = samples_rl[2] if len(samples_rl) > 2 else "NO"

    lib_dict: Dict[str, List[str]] = {}
    for lib in samples[1:]:
        lib = lib.rstrip("\n").split("\t")
        if not lib or len(lib) < 3:
            continue
        lib_dict[lib[0]] = lib[1:]

    data_dict = {}

    # reference divsum
    ref_divsum = lib_dict[ref_library][0].replace(".divsum", ".align.fam.divsum")
    nucs = int(lib_dict[ref_library][1])

    data = read_lines(Path(ref_divsum))
    s_matrix_div = data.index("-----\t------\t------\t-----------\t-------\n")
    s_matrix = data.index("Coverage for each repeat class and divergence (Kimura)\n")

    divergence: Dict[str, str] = {}
    elements_div = data[s_matrix_div + 1 : s_matrix - 2]
    for line in elements_div:
        info = line.split()
        if len(info) >= 2:
            divergence[info[1]] = info[-1]

    matrix: List[Tuple[str, List[int]]] = []
    elements = data[s_matrix + 1].split()
    for element in elements[1:]:
        element = element.split("/")[-1]
        matrix.append((element, []))

    n_el = len(matrix)
    big_row = 0.0

    # build mutable list for counts
    matrix_mut = [[name, []] for (name, _) in matrix]

    for line in data[s_matrix + 2 : s_matrix + 45]:
        info = line.split()[1:]
        if len(info) < n_el:
            continue
        for n in range(n_el):
            matrix_mut[n][1].append(int(info[n]))
        suma = sum(int(x) for x in info)
        suma_rel = suma / nucs
        if suma_rel > big_row:
            big_row = suma_rel

    family_abs = {matrix_mut[n][0]: sum(matrix_mut[n][1]) for n in range(n_el)}
    sort_abs = sorted(family_abs.items(), key=operator.itemgetter(1), reverse=True)

    # load len information from table.txt
    table_dict: Dict[str, str] = {}
    table_path = Path("table.txt")
    if table_path.exists():
        for l in read_lines(table_path):
            info = l.split()
            if len(info) >= 2:
                table_dict[info[0]] = info[1]


    # fallback monomer lengths directly from the provided monomer FASTA
    fasta_len = fasta_id_to_length(fasta_monomers)
    # ------------------------------------------------------------
    # Build a GLOBAL satellite list so ALL libraries (and subtractive
    # landscapes) share the same satellite set and ordering.
    #
    # Ordering rule:
    #   1) All satellites present in the reference divsum, ordered by
    #      decreasing abundance in the reference (original behavior).
    #   2) Satellites absent in the reference but present in any other
    #      divsum, appended and ordered by their MAX abundance across
    #      non-reference libraries.
    #
    # This satisfies the request to still *list* satellites missing in
    # the reference (they will be 0 in ref outputs), while keeping a
    # consistent order across libraries (required for rep_land
    # subtraction).
    # ------------------------------------------------------------

    # Pre-parse all non-reference libraries so we can collect "extra" satellites.
    parsed_libs: Dict[str, Dict[str, object]] = {}
    # store already-parsed reference too (for reuse later)
    parsed_libs[ref_library] = {
        "nucs": nucs,
        "divergence": divergence,
        "matrix_abs_dict": {name: vals for name, vals in matrix_mut},
        "family_abs": family_abs,
    }

    for library in lib_dict:
        if library == ref_library:
            continue
        lib_divsum = lib_dict[library][0].replace(".divsum", ".align.fam.divsum")
        nucs_lib = int(lib_dict[library][1])

        data2 = read_lines(Path(lib_divsum))
        s_matrix_div2 = data2.index("-----\t------\t------\t-----------\t-------\n")
        s_matrix2 = data2.index("Coverage for each repeat class and divergence (Kimura)\n")

        divergence2: Dict[str, str] = {}
        elements_div2 = data2[s_matrix_div2 + 1 : s_matrix2 - 2]
        for line in elements_div2:
            info = line.split()
            if len(info) >= 2:
                divergence2[info[1]] = info[-1]

        elements2 = data2[s_matrix2 + 1].split()
        matrix_mut2 = [[element.split("/")[-1], []] for element in elements2[1:]]
        n_el2 = len(matrix_mut2)

        for line in data2[s_matrix2 + 2 : s_matrix2 + 45]:
            info = line.split()[1:]
            if len(info) < n_el2:
                continue
            for n in range(n_el2):
                matrix_mut2[n][1].append(int(info[n]))
            suma = sum(int(x) for x in info)
            suma_rel = suma / nucs_lib
            if suma_rel > big_row:
                big_row = suma_rel

        family_abs2 = {matrix_mut2[n][0]: sum(matrix_mut2[n][1]) for n in range(n_el2)}
        parsed_libs[library] = {
            "nucs": nucs_lib,
            "divergence": divergence2,
            "matrix_abs_dict": {name: vals for name, vals in matrix_mut2},
            "family_abs": family_abs2,
        }

    # Reference-ordered list
    ref_ordered = [fam_name for fam_name, _ in sort_abs]
    ref_set = set(ref_ordered)

    # Collect extras + score by max abs abundance across non-ref libraries
    extra_score: Dict[str, int] = {}
    for library, pdata in parsed_libs.items():
        if library == ref_library:
            continue
        fam_abs_lib: Dict[str, int] = pdata["family_abs"]  # type: ignore[assignment]
        for fam_name, abs_count in fam_abs_lib.items():
            if fam_name in ref_set:
                continue
            prev = extra_score.get(fam_name, 0)
            if abs_count > prev:
                extra_score[fam_name] = abs_count

    extras_ordered = [
        fam for fam, _ in sorted(extra_score.items(), key=operator.itemgetter(1), reverse=True)
    ]

    all_fams_ordered = ref_ordered + extras_ordered

    defnames: List[Tuple[str, str]] = []
    eq_names: List[Tuple[str, str]] = []
    sat_digits = len(str(len(all_fams_ordered)))

    for n, fam_name in enumerate(all_fams_ordered):
        number_name = str(n + 1).zfill(sat_digits)
        monomer_len = table_dict.get(fam_name)
        if monomer_len is None:
            # satellite may be absent from reference; try monomer FASTA
            ml = fasta_len.get(fam_name)
            monomer_len = str(ml) if ml is not None else "NA"
        defnames.append((fam_name, f"{sp_name}Sat{number_name}-{monomer_len}"))
        eq_names.append((fam_name, f"{sp_name}Sat{number_name}"))

    # write equivalences.txt (now includes also satellites absent in ref)
    with Path("equivalences.txt").open("w", encoding="utf-8") as out:
        for old, new in eq_names:
            out.write(f"{old}\t{new}\n")

    family_rel_def = []
    for fam_name, defname in defnames:
        number = family_abs.get(fam_name, 0)
        rel_number = round(number / nucs, 100)
        family_rel_def.append((defname, rel_number, divergence.get(fam_name, "NA")))

    # write ref .abdiv
    with Path(f"{ref_library}.abdiv").open("w", encoding="utf-8") as out:
        for name, rel, div in family_rel_def:
            out.write(f"{name}\t{rel}\t{div}\n")

    # matrix to dict + rel
    matrix_abs_dict = parsed_libs[ref_library]["matrix_abs_dict"]  # type: ignore[assignment]
    # Determine number of divergence bins (n_div) from the reference matrix.
    # This must be defined *before* we fill missing families with zeros.
    if matrix_abs_dict:
        n_div = len(next(iter(matrix_abs_dict.values())))
    else:
        # Extremely defensive fallback: empty reference matrix (should not happen).
        n_div = 0

    matrix_rel_list = []
    for fam_name, defname in defnames:
        lista = matrix_abs_dict.get(fam_name, [0] * n_div)
        lista_rel = [x / nucs for x in lista]
        matrix_rel_list.append((defname, lista_rel))

    data_dict[ref_library] = matrix_rel_list

    # write rel matrix
    header = ["Div"] + [name for name, _ in matrix_rel_list]
    with Path(f"{ref_library}_rl.txt").open("w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for a in range(n_div):
            line = [str(a)] + [str(matrix_rel_list[b][1][a]) for b in range(len(matrix_rel_list))]
            out.write("\t".join(line) + "\n")

    # other libraries (now use the GLOBAL satellite list)
    if len(lib_dict) > 1:
        for library in lib_dict:
            if library == ref_library:
                continue

            pdata = parsed_libs[library]
            nucs_lib = int(pdata["nucs"])  # type: ignore[arg-type]
            divergence2 = pdata["divergence"]  # type: ignore[assignment]
            family_abs2 = pdata["family_abs"]  # type: ignore[assignment]
            matrix_abs_dict2 = pdata["matrix_abs_dict"]  # type: ignore[assignment]

            family_rel_def2 = []
            for fam_name, defname in defnames:
                number = family_abs2.get(fam_name, 0)
                rel_number = round(number / nucs_lib, 100)
                family_rel_def2.append((defname, rel_number, divergence2.get(fam_name, "NA")))

            with Path(f"{library}.abdiv").open("w", encoding="utf-8") as out:
                for name, rel, div in family_rel_def2:
                    out.write(f"{name}\t{rel}\t{div}\n")

            matrix_rel_list2 = []
            for fam_name, defname in defnames:
                lista = matrix_abs_dict2.get(fam_name, [0] * n_div)
                lista_rel = [round(x / nucs_lib, 100) for x in lista]
                matrix_rel_list2.append((defname, lista_rel))

            data_dict[library] = matrix_rel_list2

            with Path(f"{library}_rl.txt").open("w", encoding="utf-8") as out:
                out.write("\t".join(header) + "\n")
                for a in range(n_div):
                    line = [str(a)] + [str(matrix_rel_list2[b][1][a]) for b in range(len(matrix_rel_list2))]
                    out.write("\t".join(line) + "\n")

    # plot via R (same as original)
    for library in lib_dict:
        r_path = Path(f"{library}_rl.R")
        script = f"""library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(grid)
library(grid)
lmig <- read.table("{library}_rl.txt",header=T)
lm <- melt(lmig, id.vars=0:1)
colourCount = {len(defnames)}
ref <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)
palette1 <- rev(ref)
pdf("{library}_rl.pdf", width=11, height=7, onefile=TRUE)
ggplot(data=lm, aes(x=Div, y=value, fill=variable))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="satDNA Families", values = palette1)+
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion") +
  guides(fill=guide_legend(ncol=3, byrow=FALSE)) +
  coord_cartesian(ylim=c(0,{big_row})) +
  theme_bw() +
  theme(
    legend.position="right",
    legend.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.key.size=unit(0.35, "cm")
  )
dev.off()
"""
        write_text(r_path, script)
        try:
            run(["Rscript", str(r_path)])
        except Exception as e:
            print(f"WARNING: Could not run Rscript for {library}: {e}")

    # merged .abdiv summaries for requested subtractive pairs
    if rep_land != "NO":
        for rl in rep_land.split(","):
            rl = rl.strip()
            if not rl:
                continue
            pairs = rl.split("-")
            if len(pairs) != 2:
                continue
            merge_abdiv_pair(pairs[0], pairs[1])

    # subtractive repeat landscape
    if rep_land != "NO":
        for rl in rep_land.split(","):
            rl = rl.strip()
            if not rl:
                continue
            pairs = rl.split("-")
            if len(pairs) != 2:
                continue
            data1 = data_dict[pairs[0]]
            data2 = data_dict[pairs[1]]
            data_subs = []
            for d in range(len(data1)):
                sat1 = data1[d][1]
                sat2 = data2[d][1]
                data_subs.append((data1[d][0], [a - b for a, b in zip(sat1, sat2)]))

            with Path(f"{rl}_rl.txt").open("w", encoding="utf-8") as out:
                out.write("\t".join(header) + "\n")
                for a in range(n_div):
                    line = [str(a)] + [str(data_subs[b][1][a]) for b in range(len(data_subs))]
                    out.write("\t".join(line) + "\n")

            r_path = Path(f"{rl}_rl.R")
            script = f"""library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(grid)
subs <- read.table("{rl}_rl.txt",header=T)
s <- melt(subs, id.vars=0:1)
s1 <- subset(s,s$value>=0)
s2 <- subset(s,s$value<0)
colourCount = {len(defnames)}
ref <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)
palette1 <- rev(ref)
pdf("{rl}_rl.pdf", width=11, height=7, onefile=TRUE)
ggplot() +
  geom_bar(data=s2, aes(x=Div, y=value, fill=variable), stat="identity", position=position_stack(reverse = TRUE)) +
  geom_bar(data=s1, aes(x=Div, y=value, fill=variable), stat="identity", position=position_stack(reverse = TRUE)) +
  scale_fill_manual(name="satDNA Families", values=palette1) +
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion") +
  guides(fill=guide_legend(ncol=3, byrow=FALSE)) +
  theme_bw() +
  theme(
    legend.position="right",
    legend.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.key.size=unit(0.35, "cm")
  )
dev.off()
"""
            write_text(r_path, script)
            try:
                run(["Rscript", str(r_path)])
            except Exception as e:
                print(f"WARNING: Could not run Rscript for {rl}: {e}")





def write_original_variant_abundances(samples_file: Path, fasta_monomers: Path, out_path: Path = Path("original_variant_abundances.txt")) -> None:
    """Write one row per original FASTA sequence with abundance/divergence in primary and secondary divsum files."""
    ref_library, lib_order, parsed = parse_original_divsum_counts(samples_file)
    secondary_library = next((lib for lib in lib_order if lib != ref_library), None)
    ref_data = parsed.get(ref_library, {})
    sec_data = parsed.get(secondary_library, {}) if secondary_library else {}

    ref_rel = ref_data.get("family_rel", {}) if isinstance(ref_data, dict) else {}
    ref_div = ref_data.get("divergence", {}) if isinstance(ref_data, dict) else {}
    sec_rel = sec_data.get("family_rel", {}) if isinstance(sec_data, dict) else {}
    sec_div = sec_data.get("divergence", {}) if isinstance(sec_data, dict) else {}

    rows = []
    for rec in SeqIO.parse(str(fasta_monomers), "fasta"):
        rid = str(rec.id)
        length = len(str(rec.seq))
        rows.append((
            rid,
            length,
            float(ref_rel.get(rid, 0.0)),
            str(ref_div.get(rid, "NA")),
            float(sec_rel.get(rid, 0.0)) if secondary_library else 0.0,
            str(sec_div.get(rid, "NA")) if secondary_library else "NA",
        ))

    with out_path.open("w", encoding="utf-8") as out:
        sec_name = secondary_library or "secondary"
        out.write(f"OriginalID	Length	{ref_library}_abundance	{ref_library}_divergence	{sec_name}_abundance	{sec_name}_divergence\n")
        for rid, length, ra, rd, sa, sd in rows:
            out.write(f"{rid}	{length}	{ra}	{rd}	{sa}	{sd}\n")


def write_grouped_satellite_abundances(samples_file: Path, out_path: Path = Path("grouped_satellite_abundances.txt"), out_abs_path: Path = Path("grouped_satellite_abundances_absolute.txt")) -> None:
    """Write one row per grouped satellite after the .align.fam.divsum step, in both relative and absolute counts."""
    samples = read_lines(samples_file)
    header = samples[0].rstrip("\n").split("\t")
    ref_library = header[1]

    lib_meta: Dict[str, List[str]] = {}
    libs: List[str] = []
    for lib in samples[1:]:
        parts = lib.rstrip("\n").split("\t")
        if len(parts) >= 3:
            libs.append(parts[0])
            lib_meta[parts[0]] = parts[1:]
    secondary_library = next((lib for lib in libs if lib != ref_library), None)
    if secondary_library is None:
        print("WARNING: grouped_satellite_abundances.txt not written because no secondary library was found.")
        return

    ref_nucs = int(lib_meta[ref_library][1])
    sec_nucs = int(lib_meta[secondary_library][1])
    ref_divsum = Path(lib_meta[ref_library][0].replace('.divsum', '.align.fam.divsum'))
    sec_divsum = Path(lib_meta[secondary_library][0].replace('.divsum', '.align.fam.divsum'))

    ref_data = parse_divsum_file(ref_divsum, ref_nucs)
    sec_data = parse_divsum_file(sec_divsum, sec_nucs)
    ref_rel = ref_data["family_rel"]  # type: ignore[index]
    ref_abs = ref_data["family_abs"]  # type: ignore[index]
    ref_div = ref_data["divergence"]  # type: ignore[index]
    sec_rel = sec_data["family_rel"]  # type: ignore[index]
    sec_abs = sec_data["family_abs"]  # type: ignore[index]
    sec_div = sec_data["divergence"]  # type: ignore[index]

    eq_map: Dict[str, str] = {}
    eq_reverse: Dict[str, str] = {}
    eq_path = Path("equivalences.txt")
    if eq_path.exists():
        for line in read_lines(eq_path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                old_name, new_base = parts[0], parts[1]
                eq_map[old_name] = new_base
                eq_reverse[new_base] = old_name

    table_orig: Dict[str, str] = {}
    table_path = Path("table.txt")
    if table_path.exists():
        for line in read_lines(table_path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                table_orig[parts[0]] = parts[1]

    table_fam: Dict[str, str] = {}
    table_fam_path = Path("table.txt.fam")
    if table_fam_path.exists():
        for line in read_lines(table_fam_path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                table_fam[parts[0]] = parts[1]

    ordered = sorted(set(list(ref_rel.keys()) + list(sec_rel.keys())))

    def resolve_ids(name: str) -> Tuple[str, str]:
        # Case 1: grouped divsum still uses the original family/leader name
        if name in eq_map:
            original_id = name
            final_base = eq_map[name]
        else:
            # Case 2: grouped divsum already uses the final base or final ID with length
            candidate = name.split("#")[0]
            base = candidate.split("-")[0]
            final_base = base
            original_id = eq_reverse.get(base, "NA")
        final_len = table_fam.get(final_base) or table_orig.get(original_id)
        final_id = f"{final_base}-{final_len}" if final_len else final_base
        return original_id, final_id

    with out_path.open("w", encoding="utf-8") as out:
        out.write(f"OriginalID\tFinalID\t{ref_library}_abundance\t{ref_library}_divergence\t{secondary_library}_abundance\t{secondary_library}_divergence\n")
        for name in ordered:
            original_id, final_id = resolve_ids(name)
            out.write(f"{original_id}\t{final_id}\t{ref_rel.get(name, 0)}\t{ref_div.get(name, 'NA')}\t{sec_rel.get(name, 0)}\t{sec_div.get(name, 'NA')}\n")

    with out_abs_path.open("w", encoding="utf-8") as out_abs:
        out_abs.write(f"OriginalID\tFinalID\t{ref_library}_abundance\t{ref_library}_divergence\t{secondary_library}_abundance\t{secondary_library}_divergence\n")
        for name in ordered:
            original_id, final_id = resolve_ids(name)
            out_abs.write(f"{original_id}\t{final_id}\t{ref_abs.get(name, 0)}\t{ref_div.get(name, 'NA')}\t{sec_abs.get(name, 0)}\t{sec_div.get(name, 'NA')}\n")

def merge_abdiv_pair(lib_a: str, lib_b: str) -> None:
    """Merge two .abdiv files side by side into <lib_a>-<lib_b>.abdiv.txt."""
    def read_abdiv(path: Path):
        rows = []
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                rows.append((parts[0], parts[1], parts[2]))
        return rows

    rows_a = read_abdiv(Path(f"{lib_a}.abdiv"))
    rows_b = read_abdiv(Path(f"{lib_b}.abdiv"))

    dict_a = {name: (ab, div) for name, ab, div in rows_a}
    dict_b = {name: (ab, div) for name, ab, div in rows_b}

    ordered_names = [name for name, _, _ in rows_a]
    for name, _, _ in rows_b:
        if name not in dict_a:
            ordered_names.append(name)

    out_path = Path(f"{lib_a}-{lib_b}.abdiv.txt")
    with out_path.open("w", encoding="utf-8") as out:
        out.write(f"Satellite\t{lib_a}_abundance\t{lib_a}_divergence\t{lib_b}_abundance\t{lib_b}_divergence\n")
        for name in ordered_names:
            a_ab, a_div = dict_a.get(name, ("0", "NA"))
            b_ab, b_div = dict_b.get(name, ("0", "NA"))
            out.write(f"{name}\t{a_ab}\t{a_div}\t{b_ab}\t{b_div}\n")

def write_prefixed_copy(path: Path, prefix: str) -> None:
    """Write a copy of `path` with an extra first row whose first cell is `prefix`.

    This is meant for downstream manual inspection/merging in spreadsheets without breaking
    the original pipeline inputs/outputs.
    Output filename: <original>.prefixed
    """
    if not path.exists():
        return
    lines = read_lines(path)
    # Preserve original file exactly below the new header row
    # Use tab on the header row if the file appears tabular; otherwise single cell.
    header = prefix + ("\t" if (lines and ("\t" in lines[0])) else "")
    out_path = path.with_name(path.name + ".prefixed")
    with out_path.open("w", encoding="utf-8") as out:
        out.write(header + "\n")
        for l in lines:
            out.write(l if l.endswith("\n") else l + "\n")
# -------------------------
# main pipeline (satminer_quant)
# -------------------------

def satminer_quant(samples_file: Path, fasta_monomers: Path) -> None:
    print("Loading files.\n")

    samples = read_lines(samples_file)
    samples_rl = samples[0].rstrip("\n").split("\t")
    if len(samples_rl) < 3:
        raise ValueError(
            "SamplesFile first line must include: sp_name<TAB>ref_library<TAB>rep_land"
        )

    # Header fields
    sp_name = samples_rl[0]
    ref_library = samples_rl[1]

    lib_dict: Dict[str, List[str]] = {}
    for lib in samples[1:]:
        lib = lib.rstrip("\n")
        if not lib.strip():
            continue
        parts = lib.split("\t")
        lib_dict[parts[0]] = parts[1:]

    if ref_library not in lib_dict:
        raise KeyError(f"Reference library '{ref_library}' not found in SamplesFile.")

    ref_divsum = Path(lib_dict[ref_library][0])

    print(f"{ref_library} is the reference library.\n")
    print("Defining families.\n")

    # 1) Early table: one row per original FASTA sequence with abundance/divergence in both original divsums
    write_original_variant_abundances(samples_file, fasta_monomers)

    # original satminer_quant: divsum_ab.py ref_divsum fasta
    divsum_ab(ref_divsum, fasta_monomers)

    # Rebuild pattern.txt from the early per-sequence abundance table so secondary-only
    # variants are still represented in pattern.txt and therefore in *.align.fam.
    rebuild_pattern_from_original_abundances(Path("original_variant_abundances.txt"), fasta_monomers)

    print("Generating divsum file per families")

    # for each library, take its .divsum path -> .align, then sat_subfam2fam
    for line in samples[1:]:
        cols = line.split()
        if len(cols) < 2:
            continue
        divsum = cols[1]
        align = Path(divsum.replace(".divsum", ".align"))
        sat_subfam2fam(align, Path("pattern.txt"))

    # convert divsum to repeat landscape + plots
    divsum_to_rl(samples_file, fasta_monomers)

    # 2) Grouped table after the .align.fam.divsum step (multiple variants counted together)
    write_grouped_satellite_abundances(samples_file)

    # apply equivalences to tabular outputs
    replace_patterns(Path("table.txt"), Path("equivalences.txt"), output_suffix=".fam")

    # build final renamed FASTA outputs with the final family names and per-variant lengths
    rename_fasta_and_dim_outputs(fasta_monomers, samples_file)



def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("samples_file", nargs="?", help="Samples file (tab-separated).")
    p.add_argument("fasta_monomers", nargs="?", help="FASTA file with monomers.")
    return p


def main(argv: List[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)

    if not args.samples_file:
        print("Usage: satminer_quant_py3_allinone.py SamplesFile FastaFileMonomers")
        args.samples_file = input("Introduce Samples File: ").strip()
    if not args.fasta_monomers:
        args.fasta_monomers = input("Introduce FASTA file with monomers: ").strip()

    satminer_quant(Path(args.samples_file), Path(args.fasta_monomers))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
