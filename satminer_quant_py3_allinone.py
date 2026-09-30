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
ggplot(data=lm, aes(x=lm$Div, y=lm$value, fill=lm$variable))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name="satDNA Families", values = palette1)+
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion") +
  guides(fill=guide_legend(ncol=3, byrow=TRUE)) +
  theme(
    legend.position="right",
    legend.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.key.size=unit(0.35, "cm")
  ) +
  ylim(0,{big_row}) +
  theme_bw()
dev.off()
"""
        write_text(r_path, script)
        try:
            run(["Rscript", str(r_path)])
        except Exception as e:
            print(f"WARNING: Could not run Rscript for {library}: {e}")

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
ggplot()+
  geom_bar(data=s2,aes(x=s2$Div, y=s2$value, fill=s2$variable),stat="identity",position=position_stack(reverse = TRUE))+
  scale_fill_manual(name="satDNA Families", values=palette1)+
  labs(x="Kimura Substitution Level (%)", y="Genome Proportion")+
  theme_bw() +
  guides(fill=guide_legend(ncol=3, byrow=TRUE)) +
  theme(
    legend.position="right",
    legend.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.key.size=unit(0.35, "cm")
  ) +
  geom_bar(data=s1,aes(x=s1$Div, y=s1$value, fill=s1$variable),stat="identity",position=position_stack(reverse = TRUE))+
  scale_fill_manual(name="satDNA Families", values=palette1)+
  labs(x="Kimura Substitution Level (%)",y="Genome Proportion")+
  theme_bw()
dev.off()
"""
            write_text(r_path, script)
            try:
                run(["Rscript", str(r_path)])
            except Exception as e:
                print(f"WARNING: Could not run Rscript for {rl}: {e}")



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

    # original satminer_quant: divsum_ab.py ref_divsum fasta
    divsum_ab(ref_divsum, fasta_monomers)

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

    # apply equivalences to outputs
    replace_patterns(Path("table.txt"), Path("equivalences.txt"), output_suffix=".fam")
    replace_patterns(Path(str(fasta_monomers) + ".abc"), Path("equivalences.txt"), output_suffix=".fam")
    replace_patterns(Path(str(fasta_monomers) + ".dim.abc"), Path("equivalences.txt"), output_suffix=".fam")

    # Create separate spreadsheet-friendly copies with prefix in A1.
    # Originals are kept unchanged.
    write_prefixed_copy(Path("selection.txt"), sp_name)
    write_prefixed_copy(Path("pattern.txt"), sp_name)
    write_prefixed_copy(Path("table.txt"), sp_name)
    write_prefixed_copy(Path("equivalences.txt"), sp_name)
    write_prefixed_copy(Path("table.txt.fam"), sp_name)
    write_prefixed_copy(Path("selection.txt.extract"), sp_name)



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
