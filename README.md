# satMiner-quant

A Python 3 implementation of the `satminer_quant` workflow for quantifying and analysing satellite DNA repeat landscapes.

The pipeline integrates the Python helper scripts used by the original workflow into a single executable and produces family-level alignments, divergence files, repeat landscapes, abundance/divergence tables, and consistently renamed FASTA files.

## Features

- Python 3 implementation of the original workflow.
- Defines satellite families from monomer FASTA sequences.
- Retains variants that are exclusive to the secondary library.
- Calculates abundance and divergence for individual variants before family-level processing.
- Rebuilds family-level alignments and `divsum` files.
- Generates repeat-landscape PDF plots.
- Produces relative and absolute abundance tables.
- Generates mappings between original and final sequence IDs.
- Renames both monomer FASTA and DIM FASTA outputs consistently.
- Uses monomer length in final satellite IDs, including for DIM FASTA output.
- Keeps plot dimensions and legend layout consistent.

## Requirements

The workflow requires:

- Python 3
- Biopython
- R / `Rscript`
- Perl
- `calcDivergenceFromAlign.pl`
- The R packages required by the generated plotting scripts, including `ggplot2` and the packages used by the original workflow.

Install Biopython with:

```bash
pip install biopython
```

`calcDivergenceFromAlign.pl` must be available in `$PATH`.

Check the main external dependencies with:

```bash
python3 --version
Rscript --version
perl --version
calcDivergenceFromAlign.pl
```

## Usage

Make the script executable:

```bash
chmod +x satminer_quant_py3_allinone_v29.py
```

Run:

```bash
./satminer_quant_py3_allinone_v29.py SAMPLES_FILE MONOMER_FASTA
```

Example:

```bash
./satminer_quant_py3_allinone_v29.py samples_dvit.txt dvit_sat.fasta
```

or:

```bash
python3 satminer_quant_py3_allinone_v29.py samples_dvit.txt dvit_sat.fasta
```

## Input files

### Samples file

The samples file defines the libraries and the repeat-landscape comparisons.

The first line contains:

```text
sp_name    ref_library    rep_land
```

where:

- `sp_name` is the species/sample prefix.
- `ref_library` is the primary/reference library.
- `rep_land` defines the repeat-landscape comparison(s).

The following lines contain the library-specific information required by the workflow.

### Monomer FASTA

The second input is a FASTA file containing the satellite monomer sequences.

Example:

```text
>dvitR1CL100
ACGT...
>dvitR2CL200
ACGT...
```

The original FASTA IDs are retained in the intermediate mapping files and are subsequently converted to final satellite IDs.

---

# Protocol

## 1. Load the input data

The pipeline reads:

- the samples file;
- the reference/primary library;
- the secondary library or libraries;
- the monomer FASTA.

The reference library is reported at the beginning of the run.

## 2. Generate the initial `divsum` data

The initial alignments/divsum data are generated before the family-level modifications.

This stage provides the abundance and divergence information for the individual sequences in the original FASTA.

## 3. Generate `original_variant_abundances.txt`

Before modifying the alignments, the pipeline creates:

```text
original_variant_abundances.txt
```

This contains one row for every sequence in the original FASTA.

The table includes:

```text
OriginalID
Length
<primary>_abundance
<primary>_divergence
<secondary>_abundance
<secondary>_divergence
```

This is an important intermediate file because a sequence may be absent from the primary divsum but present in the secondary library.

## 4. Define families and variants

Family relationships are then constructed using the original sequence information.

The generated:

```text
pattern.txt
```

contains the relationships used to modify the alignments.

Variants exclusive to the secondary library are retained.

For example, if:

```text
F1V2dvR5CL397b_A1
F1V1dvR5CL397a_A1
```

belong to the same original variant group, both can be included in `pattern.txt`, even if the second sequence has zero abundance in the primary library.

This allows the secondary-specific variant to enter the subsequent family-level alignment.

## 5. Modify alignments

The family definitions are applied to the relevant alignment files.

Family-level alignments are generated as:

```text
*.align.fam
```

The variants belonging to a satellite family are represented together at this stage.

## 6. Recalculate divergence

For every family-level alignment, the pipeline runs:

```text
calcDivergenceFromAlign.pl
```

to create:

```text
*.align.fam.divsum
```

These are the divsum files used for the final repeat landscapes and grouped satellite measurements.

## 7. Generate repeat landscapes

R scripts are generated and executed for the individual libraries and for the comparisons defined in `rep_land`.

Typical outputs:

```text
<library>_rl.R
<library>_rl.pdf

<library1>-<library2>_rl.R
<library1>-<library2>_rl.pdf
```

The plots use a fixed PDF size and a consistent legend arrangement.

The legend is arranged by columns rather than by rows.

## 8. Generate grouped abundance/divergence tables

After family-level divsum generation, the pipeline creates:

```text
grouped_satellite_abundances.txt
grouped_satellite_abundances_absolute.txt
```

These tables contain one row per final satellite.

The columns include:

```text
OriginalID
FinalID
<primary>_abundance
<primary>_divergence
<secondary>_abundance
<secondary>_divergence
```

The rows are sorted alphabetically by `FinalID`.

The `_absolute` file contains absolute read counts instead of relative abundance.

## 9. Generate final variant tables

The pipeline also produces:

```text
final_variant_abundances.txt
final_variant_abundances_absolute.txt
```

These retain the individual sequence level.

The tables contain:

```text
OriginalID
FinalID
<primary>_abundance
<primary>_divergence
<secondary>_abundance
<secondary>_divergence
```

The absolute version reports counts rather than relative abundance.

## 10. Generate equivalence mappings

The pipeline generates:

```text
equivalences.txt
equivalences.txt.fam
```

`equivalences.txt.fam` documents the mapping between the original FASTA sequence IDs and their final IDs.

For example:

```text
F14V3dvR1CL81_B1    DviSat10A-211
F14V4dvFR2CL249_B1  DviSat10B-211
```

This file is useful for tracing every renamed sequence back to its original identifier.

## 11. Rename the FASTA outputs

The final monomer FASTA is written as:

```text
*.fasta.fam
```

The corresponding DIM FASTA is:

```text
*.fasta.dim.fam
```

Final IDs have the form:

```text
DviSat10A-211#Satellite/DviSat10A-211
```

The length encoded in the final ID is the **monomer length**, not necessarily the actual length of the DIM sequence.

Therefore, the length component is consistent between:

```text
*.fasta.fam
*.fasta.dim.fam
```

Variants in the same final satellite can receive suffixes such as:

```text
DviSat01A-100
DviSat01B-105
DviSat01C-98
```

A sequence exclusive to the secondary library can receive its own satellite ID rather than incorrectly inheriting a suffix from a satellite shared with the primary library.

---

# Main output files

A typical analysis produces some or all of:

```text
original_variant_abundances.txt
final_variant_abundances.txt
final_variant_abundances_absolute.txt
grouped_satellite_abundances.txt
grouped_satellite_abundances_absolute.txt

equivalences.txt
equivalences.txt.fam

pattern.txt
selection.txt
table.txt
table.txt.fam
selection.txt.extract

*.align.fam
*.align.fam.divsum
*.abdiv

*_rl.R
*_rl.pdf

*.fasta.fam
*.fasta.dim.fam
```

The exact files depend on the number of libraries and comparisons specified in the samples file.

## Important intermediate files

| File | Purpose |
|---|---|
| `original_variant_abundances.txt` | Abundance/divergence of every original FASTA sequence before family-level modification |
| `pattern.txt` | Family/variant relationships used to modify alignments |
| `equivalences.txt` | Original family/sequence naming relationships |
| `equivalences.txt.fam` | Mapping from original sequence IDs to final IDs |
| `final_variant_abundances.txt` | Final names plus relative abundance/divergence per sequence |
| `final_variant_abundances_absolute.txt` | Same information using absolute counts |
| `grouped_satellite_abundances.txt` | Final satellite/family-level relative abundance/divergence |
| `grouped_satellite_abundances_absolute.txt` | Final satellite/family-level absolute counts |
| `*.align.fam` | Family-level alignments |
| `*.align.fam.divsum` | Divergence results after family-level grouping |
| `*.fasta.fam` | Final renamed monomer FASTA |
| `*.fasta.dim.fam` | Final renamed DIM FASTA |

---

# Example workflow

```bash
./satminer_quant_py3_allinone_v29.py samples_dvit.txt dvit_sat.fasta
```

Expected high-level progress:

```text
Starting satMiner quantification...

  - Loading samples file
  - Loading monomer FASTA

Reference library: dvit_female

Step 1/4 - Defining families...

Step 2/4 - Generating family-level divsum files...

Step 3/4 - Building repeat landscapes and summary tables...

Step 4/4 - Renaming final FASTA outputs...

Done.
```

External programs are run quietly during successful execution. Their error output is displayed when a command fails.

---

# Handling secondary-specific variants

A key feature is retaining sequences found only in the secondary library.

For example:

```text
F1V2dvR5CL397b_A1
F1V1dvR5CL397a_A1
```

may have:

```text
F1V2dvR5CL397b_A1    primary > 0    secondary > 0
F1V1dvR5CL397a_A1    primary = 0    secondary > 0
```

The second sequence must still be included in the family definition and alignment processing.

The early `original_variant_abundances.txt` table is used to retain this information before the family-level divsum files are regenerated.

---

# Reproducibility

For reproducible analyses, keep together:

```text
satminer_quant_py3_allinone_v29.py
samples_*.txt
*_sat.fasta
calcDivergenceFromAlign.pl
```

and record software versions:

```bash
python3 --version
Rscript --version
perl --version
```

Also record the versions of the relevant Python and R packages.

---

# Relationship to the original workflow

This project is a Python 3 modernization and integration of the workflow originally implemented using Python 2 helper scripts.

The pipeline integrates the required Python components so that the main analysis can be run from a single Python 3 executable, while external Perl and R tools remain required for divergence calculation and plotting.

If this repository incorporates code from an upstream project, retain the upstream copyright notices and license terms.

