# Contributing

Contributions and reproducible bug reports are welcome after the maintainers publish the repository's licence and contribution policy.

## Before changing code

1. Create a focused branch and keep generated biological results out of the commit.
2. Describe the input type, operating system, Python/R versions, package versions, and external command versions relevant to the change.
3. Do not commit patient data, sensitive sample metadata, unpublished assemblies, credentials, or large generated databases.

## Local checks

From an activated environment:

```bash
python -m py_compile *.py
python PhyloCode.py --help
python cgMLSTtoPhyloCode.py --help
```

For workflow changes, also run a small, redistributable fixture through both entry points and compare:

- retained/excluded allele counts;
- the per-locus processing summary and log;
- generated JSON keys and numerical values;
- converted CSV shape, missing-value count, and distance matrix;
- map generation, when plotting code changed.

There is currently no automated test suite. New bug fixes should preferably include a minimal test that fails before the fix and passes afterwards.

## Style and scope

- Preserve input and output compatibility unless a breaking change is explicitly documented.
- Use `subprocess.run` with an argument list and `check=True` for external commands.
- Validate paths, matrices, identifiers, and empty inputs before long-running work.
- Add user-facing options to `--help`, README examples, and the input/output reference.
- Avoid drive-by reformatting of unrelated files.

## Bug reports

Include the exact command (with sensitive paths and names redacted), the traceback or relevant log section, and a minimal synthetic input if possible. State whether the run was fresh or resumed from an existing locus directory.

