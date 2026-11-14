# E3 — Implementation

The frozen scientific plan is [E3.md](E3.md). This directory contains the experiment code, **not a generated research report**.

| File | What it does |
|---|---|
| `config.yaml` | Defines all 84 covariance configurations and seeds |
| `design.py` | Creates population spectra, experiment grid, random streams |
| `simulation.py` | Samples covariance matrices and measures losses/predictions |
| `analysis.py` | Computes numerical estimates, intervals and secondary metrics |
| `figures.py` | Plots the seven diagnostic PNG figures |
| `run.py` | Provides safe commands for testing, verification and new runs |

Reusable algorithms are in [src/spectral_stability](../../src/spectral_stability/), with tests in [tests](../../tests/). Saved evidence is in [results/e3](../../results/e3/).

## Check the original pilot

From the repository root, with Python 3.12 and `requirements-e3.txt` installed:

```sh
python experiments/e3/run.py validate
python experiments/e3/run.py verify
```

- `validate` runs automated tests **without modifying** tracked records.
- `verify` checks the archived pilot's immutable data, hashes and numerical summaries **without modifying or resimulating** anything.
- Neither command writes a manuscript.

## Reproduce into a new directory

```sh
python experiments/e3/run.py smoke --output results/e3/reproductions/smoke-1
python experiments/e3/run.py pilot --output results/e3/reproductions/pilot-1
```

The commands refuse existing output directories. Both run tests first; `pilot` also runs four smoke reproducibility checks before the 84 configurations. Numerical summaries and **PNG-only** figures are generated automatically; written reports are not.

Use `python experiments/e3/run.py analyze --output PATH` to regenerate results **only for a newly created run with the current numerical source fingerprints**. The archived pilot is intentionally read-only.

The independent 5,000-replication confirmatory stage has not been implemented or run. The old run's source and software identifiers remain in its original manifest and Git history; they are historical provenance, not hashes of today's cleaned code.
