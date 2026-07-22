# Inputs and outputs

This document describes behaviour visible in the current scripts. It is not a biological validation protocol.

## `PhyloCode.py`

### Allele FASTA (`--in`)

- The file base name is used as the locus name. For `abcZ.fasta`, the locus is `abcZ`.
- At least three sequences are required before and after BLAST filtering.
- Each identifier must be `<locus>_<allele_number>`, where the final component is an integer.
- Identifiers must be unique and sequences must not be empty.
- The BLAST filtering step retains hits with the requested identity threshold and approximately full query coverage. Consult the log and `*_excluded_blast.tsv` before accepting a database.

### Reference BLAST database (`--db`)

Pass the database prefix produced by `makeblastdb`, not merely the source FASTA path. It should represent complete genomes of the intended target species because it is used to exclude alleles outside the target group.

### Options

| Option | Meaning | Default |
|---|---|---:|
| `-i`, `--in` | Input multi-FASTA | required |
| `-d`, `--db` | Nucleotide BLAST database prefix | required |
| `-o`, `--out` | Output root | current directory |
| `--identity` | Minimum percent identity used in filtering | `95` |
| `--force` | Remove and rebuild an existing locus work directory | off |
| `--discard_temp` | Delete intermediate files after the run | off |
| `--method` | `FastTree` or `RAxML` | `FastTree` |
| `--threads` | RAxML-NG thread count | `4` |
| `--max-lpnet-taxa` | Maximum medoids passed to lpnet (1-80) | `80` |

### Output layout

```text
<output>/
|-- DATABASE/
|   `-- <locus>.json
`-- <locus>/
    |-- <locus>.log
    |-- <locus>.eps
    |-- <locus>_represent.fasta
    |-- <locus>_processing_summary.tsv
    `-- intermediate alignments, matrices, trees, lists and Nexus files
```

`DATABASE/<locus>.json` maps allele identifiers such as `abcZ_17` to numerical angular values. `cgMLSTtoPhyloCode.py` consumes all JSON files in that directory.

The locus log is also used for resumption: stages logged as `End` may be skipped on a later invocation. Use `--force` when a clean rebuild is required.

## `cgMLSTtoPhyloCode.py`

### Locus-name table (`--locus`)

This is a tab-delimited, headerless text file:

```text
valid_locus_name	optional_input_name
```

The first column determines output order and must match JSON file names in `--db`. If the input profile already uses that name, the second column can be omitted.

### pubMLST input (`--pubmlst` or `--list`)

`--pubmlst` accepts one Excel workbook. `--list` accepts a plain-text file containing one Excel workbook path per line. For each workbook:

- the first two columns are read;
- the first row is treated as a header;
- the columns represent locus label and allele number;
- `-` is treated as missing (`-1`);
- the sample name is the workbook file name without its extension.

Despite the legacy CLI help text mentioning CSV in one place, the current implementation uses `pandas.read_excel`.

### EnteroBase input (`--enterobase`)

The current parser expects a headerless, tab-delimited matrix in which:

- row 1, columns 3 onward contain locus names;
- subsequent rows contain the sample name in column 1;
- columns 3 onward contain allele numbers;
- `-` is treated as missing (`-1`).

Confirm an export against a small expected result before processing a large study.

### Previously converted input (`--phylocode`)

The script reads a CSV with locus names in the first column and samples in subsequent columns. This mode is for redrawing or re-referencing an existing PhyloCode table and does not require `--db` or `--locus` unless their information is needed for the selected display.

### Missing and unknown values

- Missing input alleles become `-1`.
- An allele number absent from a locus JSON initially becomes an unknown value and is ultimately written as `-1` unless optional nearest-allele recovery succeeds.
- In ordinary maps, missing values are displayed as black.
- Distance calculation replaces missing `-1` values with `-3.14` before summing absolute per-locus differences. This behaviour should be considered when interpreting distances between profiles with different amounts of missing data.

### Optional nearest-allele recovery

For unknown alleles, `--genome` and `--represent` can search a sample assembly against representative alleles:

- `--genome` is a whitespace-delimited file containing `sample_name genome.fasta` per line.
- `--represent` is a directory containing `<locus>_represent.fasta` files.
- BLAST hits require 99% identity and 99% query coverage in the current implementation.
- Details are written to `<prefix>_alternate_allele.xlsx` when applicable.

### Principal outputs

| Output | Description |
|---|---|
| `<prefix>_PhyloCode.csv` | Converted locus-by-sample values |
| `<prefix>_distance_matrix.csv` | Sum of absolute per-locus differences between samples |
| `<prefix>_alternate_allele.xlsx` | Optional nearest-allele BLAST details |
| `<prefix>_map.eps` | Optional linear or circular colour map |

Use `--map` for a linear map. Use `--circular N` to request a circular map when the number of strains is at most `N`, falling back to linear above that threshold. `--ref SAMPLE` expresses values relative to a reference sample; `--remainref` keeps that sample in the plot. `--strain FILE` limits the displayed strains, and `--label_interval` controls locus label spacing.

