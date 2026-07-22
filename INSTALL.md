# Installation

PhyloCode combines Python, R, and native bioinformatics programs. A Linux environment is recommended because that is the environment assumed by the command names in the scripts.

## 1. Install system software

On Debian or Ubuntu, the following is a starting point (package availability and names vary by release):

```bash
sudo apt update
sudo apt install \
  python3 python3-venv python3-dev build-essential \
  r-base r-base-dev libglpk-dev \
  ncbi-blast+ mafft snp-dists fasttree raxml-ng \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev
```

RAxML-NG is only needed when `--method RAxML` is selected. Evince is only needed for `--open_viewer`.

## 2. Create the Python environment

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`rpy2` must be built against, and run with, the same R installation. If R is installed in a non-standard location, configure `R_HOME` before installing `requirements.txt`.

## 3. Install the R packages

Start R in the same environment and run:

```r
install.packages(c("Rglpk", "devtools"))
devtools::install_github("yukimayuli-gmz/lpnet", ref = "main")
```

`Rglpk` requires the GLPK development library when installed from source. The upstream [`lpnet` repository](https://github.com/yukimayuli-gmz/lpnet) documents its GitHub installation; [`Rglpk`](https://cran.r-project.org/package=Rglpk) documents the GLPK system requirement.

## 4. Verify the environment

```bash
python --version
R --version
blastn -version
makeblastdb -version
mafft --version
snp-dists --version
FastTree 2>&1 | head -n 1
raxml-ng --version

python -c "import Bio, jinja2, matplotlib, numpy, openpyxl, pandas, rpy2, scipy"
python -c "import rpy2.robjects.packages as p; p.importr('lpnet'); p.importr('Rglpk')"
python PhyloCode.py --help
python cgMLSTtoPhyloCode.py --help
```

If only FastTree is used, a missing `raxml-ng` is acceptable. Conversely, `FastTree` is not required when all runs explicitly use `--method RAxML`.

## Platform notes

- The scripts call `FastTree` with that exact capitalization.
- `--open_viewer` launches `evince`; omit it on systems without Evince and open the EPS output manually.
- Windows may require extra compiler and R configuration for `rpy2`, `Rglpk`, and `lpnet`. WSL2 or a Linux container is likely simpler, but has not been validated in this repository.
- Version pins are intentionally not claimed because a known-good environment has not yet been recorded. Before a release, capture tested versions and create a lock file or environment specification.

