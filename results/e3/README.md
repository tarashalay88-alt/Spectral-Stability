# E3 — Exploratory pilot evidence

**Status:** 84 configurations × 1,000 replications; independent confirmation pending. This folder stores data and figures, **not a scientific report**.

| Location | What it contains |
|---|---|
| [pilot/figures/](pilot/figures/) | Seven diagnostic PNG figures |
| [pilot/summary.csv](pilot/summary.csv) | Primary numerical results for the 84 configurations |
| [pilot/contrasts.csv](pilot/contrasts.csv) | Two-background controlled comparisons |
| [pilot/ranking.csv](pilot/ranking.csv) | Secondary risk ordering |
| [pilot/calibration.csv](pilot/calibration.csv) | Secondary capped-estimate calibration |
| [pilot/replications.csv](pilot/replications.csv) | All 84,000 original matched observations |
| [pilot/run_manifest.json](pilot/run_manifest.json) | Original settings, seeds, environment and source hashes |

The original `validation.json` and `validation.log` record the 32-test gate; `smoke/replications.csv` and `smoke/run_manifest.json` record the separate 160-row engineering smoke test, **not another research experiment**.

The historical pilot manifest lists optional SVG images and automatically generated prose that have been removed from today's repository. Those historical outputs remain available in [the original archival commit](https://github.com/tarashalay88-alt/Spectral-Stability/tree/24c3f914ca0b74ba914b7bd47cf26e201707eb5b/results/e3). No raw observation, primary numerical summary, retained PNG, or frozen protocol was modified.

To verify the *currently retained* evidence without generating files, run `python experiments/e3/run.py verify` from the repository root. It independently recalculates the main numerical summaries and checks original hashes; it does not claim to rerun the historical environment.

See the [frozen protocol](../../experiments/e3/E3.md) and [implementation guide](../../experiments/e3/README.md).
