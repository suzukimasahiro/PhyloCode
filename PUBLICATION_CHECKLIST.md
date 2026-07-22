# GitHub publication checklist

The code and basic user documentation are present, but the following owner decisions are still required before making the repository public.

## Required decisions

- [ ] Choose a repository name and add the final GitHub URL to the README.
- [x] Add an MIT `LICENSE` file. Before publication, confirm that all copyright holders and the institution permit this licence.
- [ ] Add author names, affiliations, ORCIDs, preferred citation, repository URL, and release DOI to `CITATION.cff`.
- [ ] State the supported operating systems and the exact tested versions of Python, R, packages, and external tools.
- [ ] Confirm the intended scientific interpretation of PhyloCode values and distance calculations with the method authors.
- [ ] Decide whether the clinical/research-use disclaimer in the README is sufficient under institutional policy.

## Reproducibility and quality

- [ ] Add a small synthetic or openly licensed example dataset with expected outputs.
- [ ] Add automated tests for FASTA validation, medoid selection, profile parsing, conversion, and missing-value distance behaviour.
- [ ] Add continuous integration once a reproducible environment is available.
- [ ] Run the complete FastTree workflow on the example data.
- [ ] Run the complete RAxML-NG workflow, or clearly mark it as untested.
- [ ] Verify that `--discard_temp`, resume behaviour, reference adjustment, circular plotting, and nearest-allele recovery work as documented.
- [ ] Pin or lock a known-good environment for a release.

## Data and security

- [ ] Remove human/patient identifiers, private sample metadata, credentials, absolute local paths, and unpublished sequence data.
- [ ] Confirm that all included example sequences and third-party assets can legally be redistributed.
- [ ] Review the full Git history for material that should not become public.
- [ ] Enable GitHub secret scanning and dependency alerts where appropriate.

## Repository presentation

- [ ] Add a concise repository description and topics such as `cgmlst`, `phylogenetics`, and `bioinformatics`.
- [ ] Add a release/changelog policy and create the first tagged release only after an end-to-end validated run.
- [ ] Configure issue templates for bug reports and feature requests if external contributions will be accepted.
- [ ] Add a support/contact route and a code of conduct if the project will accept community participation.

## Suggested final checks

```bash
git status --short
git diff --check
python -m py_compile *.py
python PhyloCode.py --help
python cgMLSTtoPhyloCode.py --help
```

The directory is initialized as a Git work tree on the `main` branch. Before the first commit, confirm that the untracked-file list contains no local, sensitive, or generated data, then configure the intended author identity and GitHub remote.
