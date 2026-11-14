# E3 — Simulation results archive

**Stage:** Exploratory pilot, run on 2026-11-13. The independent 5,000-replication confirmation has **not** been run.

The frozen scientific protocol is [experiments/e3/E3.md](../../experiments/e3/E3.md). This directory archives the original pilot evidence and the verification materials. It is separate from the eventual human-written results and discussion.

## Read first

- [pilot/summary.csv](pilot/summary.csv) — one row for each of the **84** covariance configurations; estimated risk, oracle and plug-in errors, confidence intervals, and reliability labels
- [pilot/figures/](pilot/figures/) — seven diagnostic plots supplied in PNG and SVG
- [pilot/pilot_report.md](pilot/pilot_report.md) — Codex's **technical** exploratory summary (not the final research report)
- [pilot/contrasts.csv](pilot/contrasts.csv) — 42 background-spectrum comparisons
- [pilot/ranking.csv](pilot/ranking.csv) and [pilot/calibration.csv](pilot/calibration.csv) — secondary analyses

## Complete provenance and reproducibility

- [pilot/replications.csv](pilot/replications.csv) — all **84,000** matched replication records, including loss, raw and capped plug-in prediction, and sample boundary gap
- [pilot/run_manifest.json](pilot/run_manifest.json) — configuration, seeds, exact source and artifact hashes, software versions, and run metadata
- [pilot/audit.json](pilot/audit.json), [pilot/postflight_report.md](pilot/postflight_report.md), [pilot/regeneration_reference.json](pilot/regeneration_reference.json), and [postflight_audit.py](postflight_audit.py) — post-run checks and regeneration evidence
- [validation.json](validation.json), [validation.log](validation.log), and [smoke/](smoke/) — prerequisite tests and smoke-run evidence
- [pilot/figure_index.json](pilot/figure_index.json) — links saved figures to the run

The implementation is in [experiments/e3](../../experiments/e3/), mathematical code in [src/spectral_stability](../../src/spectral_stability/), and tests in [tests](../../tests/). Source provenance: Codex commit `5ecd5c9b6930ec2b2821483e8855bb9aeb7255c5` and run manifest's source fingerprint.

## Interpretation policy

These are **exploratory** observations from 1,000 replications per configuration, not confirmed findings or universal reliability boundaries. The frozen protocol has not been amended. Review the pilot and the analysis methods before independently running confirmation. Do not infer that 32 passing automated tests or a successful postflight audit proves statistical validity.

Only the canonical source and results trees were imported; the duplicate `deliverables/e3-pilot/` ZIP, copies of figures, and patch were intentionally not duplicated in `main`.
