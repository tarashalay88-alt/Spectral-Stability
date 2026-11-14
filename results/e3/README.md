# E3 — Pilot data and figures

**Stage:** Exploratory pilot only (84 configurations × 1,000 replications). Independent confirmation has not been run.

## Explore the results

- [`pilot/summary.csv`](pilot/summary.csv) — one row per configuration: risk estimates, oracle and plug-in errors, confidence intervals, and reliability labels.
- [`pilot/figures/`](pilot/figures/) — seven pilot figures, in **PNG only**.
- [`pilot/contrasts.csv`](pilot/contrasts.csv) — background-spectrum comparisons at fixed sample size and boundary gap.
- [`pilot/ranking.csv`](pilot/ranking.csv) and [`pilot/calibration.csv`](pilot/calibration.csv) — secondary numerical analyses.

## Reproduce and verify

- [`pilot/replications.csv`](pilot/replications.csv) — all **84,000** original matched simulation records (not modified).
- [`pilot/run_manifest.json`](pilot/run_manifest.json) — the original simulation's parameters, random seeds, source hashes and software versions.
- [`pilot/audit.json`](pilot/audit.json), [`pilot/regeneration_reference.json`](pilot/regeneration_reference.json), [`validation.json`](validation.json), and [`smoke/`](smoke/) — recorded validation and smoke-test evidence.
- [`pilot/figure_index.json`](pilot/figure_index.json) — original output index, including hashes for historical files.

The frozen protocol is [`experiments/e3/E3.md`](../../experiments/e3/E3.md); implementation and tests are in [`experiments/e3/`](../../experiments/e3/), [`src/spectral_stability/`](../../src/spectral_stability/), and [`tests/`](../../tests/).

**Original-run integrity:** The manifests and audit describe the **complete original run**, including automatically generated reports and SVG figures. Those presentation duplicates are omitted from the current repository tree at the author's request. They remain recoverable from [the original E3 import commit](https://github.com/tarashalay88-alt/Spectral-Stability/tree/24c3f914ca0b74ba914b7bd47cf26e201707eb5b/results/e3). As a result, a full-artifact audit of this curated tree alone will not pass without restoring the omitted outputs. The simulation records, metrics, PNG figures and frozen protocol have not been changed.

No author-written findings or final report have been added. Results are preliminary and must be reviewed before confirmation.
