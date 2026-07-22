# PhyloCode

PhyloCode is a command-line research workflow for converting allele sequences and cgMLST profiles into circular, locus-specific numerical codes. It contains two user-facing programs:

- `PhyloCode.py` builds a JSON conversion database for one locus from an allele FASTA file.
- `cgMLSTtoPhyloCode.py` applies those JSON databases to pubMLST or EnteroBase cgMLST profiles, calculates pairwise distances, and can draw linear or circular colour maps.

> [!IMPORTANT]
> This repository currently contains research software. Its biological assumptions and outputs should be independently validated before clinical, diagnostic, or public-health use.

## Workflow

```text
allele FASTA + species BLAST database
                  |
                  v
             PhyloCode.py
                  |
                  v
       DATABASE/<locus>.json
                  |
                  + cgMLST profile + locus-name table
                  |
                  v
       cgMLSTtoPhyloCode.py
                  |
                  +-- <prefix>_PhyloCode.csv
                  +-- <prefix>_distance_matrix.csv
                  +-- <prefix>_map.eps (optional)
```

## Requirements

Linux is the recommended environment. The workflow requires Python 3, R, and several command-line bioinformatics tools.

- Python packages: Biopython, Jinja2 3.x, matplotlib, NumPy, openpyxl, pandas, rpy2, and SciPy
- Commands: NCBI BLAST+ (`blastn`, `makeblastdb`), MAFFT, `snp-dists`, and either FastTree or RAxML-NG
- R packages: `lpnet` and `Rglpk`
- Optional: Evince, only for `cgMLSTtoPhyloCode.py --open_viewer`

See [INSTALL.md](INSTALL.md) for installation and environment checks.

## Quick start: build a locus database

The input FASTA file name (without its extension) becomes the locus name. It must contain at least three non-empty sequences. Every identifier must be unique and use the form `<locus>_<integer>`.

For example, `abcZ.fasta` should look like:

```fasta
>abcZ_1
ATGC...
>abcZ_2
ATGT...
>abcZ_3
ATGA...
```

Create a nucleotide BLAST database containing complete reference genomes for the target species:

```bash
makeblastdb -in target_species_genomes.fasta -dbtype nucl -out reference/target_species
```

Run PhyloCode:

```bash
python PhyloCode.py \
  --in abcZ.fasta \
  --db reference/target_species \
  --out results \
  --identity 95 \
  --method FastTree
```

The main reusable outputs are:

- `results/DATABASE/abcZ.json`: allele-to-PhyloCode conversion map
- `results/abcZ/abcZ_represent.fasta`: representative allele sequences
- `results/abcZ/abcZ.eps`: locus map
- `results/abcZ/abcZ_processing_summary.tsv`: processing summary
- `results/abcZ/abcZ.log`: resumable processing log

Run the command once per locus while using the same output directory to accumulate JSON files under `DATABASE/`. Use `--force` to replace an existing locus work directory and `--discard_temp` to remove intermediate files after success.

To use RAxML-NG instead of FastTree:

```bash
python PhyloCode.py -i abcZ.fasta -d reference/target_species \
  -o results --method RAxML --threads 8
```

Show all options with `python PhyloCode.py --help`.

## Quick start: convert cgMLST profiles

Prepare a tab-delimited locus-name table. The first column is the database/valid locus name; the optional second column is an alternative name used by the input profile:

```text
abcZ
adk	ADK_ALT_NAME
aroE
```

For pubMLST input, supply either one Excel workbook with `--pubmlst` or a text file containing one workbook path per line with `--list`. The first two workbook columns are read as locus and allele number.

```bash
python cgMLSTtoPhyloCode.py \
  --pubmlst sample.xlsx \
  --locus loci.tsv \
  --db results/DATABASE \
  --out comparison \
  --map
```

For an EnteroBase profile exported in the layout expected by the script:

```bash
python cgMLSTtoPhyloCode.py \
  --enterobase enterobase_profile.tsv \
  --locus loci.tsv \
  --db results/DATABASE \
  --out comparison \
  --circular 20
```

To redraw a map from an already converted file:

```bash
python cgMLSTtoPhyloCode.py \
  --phylocode comparison_PhyloCode.csv \
  --out comparison_redrawn \
  --map
```

See [docs/INPUTS_AND_OUTPUTS.md](docs/INPUTS_AND_OUTPUTS.md) for the accepted formats, missing-value behaviour, optional nearest-allele search, and output definitions.

## Reproducibility notes

- Record the versions of Python, R, all Python/R packages, and external commands used for a run.
- The selected phylogeny method affects the conversion database. FastTree is seeded with `1234`; RAxML-NG is seeded with `1` and uses 100 bootstrap trees.
- Keep the complete per-locus directory and log until the database has been checked. `--discard_temp` is intended for a completed, verified run.
- Input and output paths should not contain ambiguous shell metacharacters. Invoke the scripts from the repository root so their sibling modules can be imported.

## Development

The codebase currently consists of standalone Python scripts rather than an installable Python package. See [CONTRIBUTING.md](CONTRIBUTING.md) for validation and contribution guidance.

## Citation and licence

This software is distributed under the [MIT License](LICENSE). Citation metadata has not yet been supplied. Before making the repository public, the maintainers should add author, ORCID, publication, and repository details in `CITATION.cff`. See [PUBLICATION_CHECKLIST.md](PUBLICATION_CHECKLIST.md).
