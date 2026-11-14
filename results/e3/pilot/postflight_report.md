# E3 completion and postflight report

Implementation and frozen exploratory pilot completed on 2026-11-13.
The authoritative `experiments/e3/E3.md` is byte-for-byte unchanged from
`afa3dcbdb54909d82193c4a9d6f78d785d5b62db`. No confirmation was run. No Git commit or remote push
was made; the source fingerprint is `e5b01b6e772e363418487848636f9861b78bd25724304b911e7164c67c26d1c9`.

## Files created

- `.gitignore`
- `requirements-e3.txt`
- `experiments/e3/README.md`
- `experiments/e3/analysis.py`
- `experiments/e3/config.yaml`
- `experiments/e3/design.py`
- `experiments/e3/figures.py`
- `experiments/e3/report.py`
- `experiments/e3/run.py`
- `experiments/e3/simulation.py`
- `src/spectral_stability/__init__.py`
- `src/spectral_stability/e3_math.py`
- `src/spectral_stability/e3_statistics.py`
- `tests/__init__.py`
- `tests/test_e3_math.py`
- `tests/test_e3_statistics.py`

Generated outputs are in `results/e3/`: `validation.json`, `validation.log`,
`smoke/` and `pilot/`. Each run has `replications.csv`, `summary.csv`,
`contrasts.csv`, `ranking.csv`, `calibration.csv`, `run_manifest.json`,
`figure_index.json`, `pilot_report.md`, and seven figures in both PNG and SVG.
The pilot additionally has this postflight report, `audit.json`, and
`regeneration_reference.json`; `results/e3/postflight_audit.py` reproduces the
postflight audit. Unpopulated secondary smoke tables represent incomplete
seven-gap strata, not missing pilot simulations.

## Implementation and validation

Existing Gaussian sampling, centered covariance, top-k basis and oracle helpers
are reused. E3 records unclipped Frobenius projector loss and unmodified raw
cross-eigenvalue sums with denominator k(n−1). Capping is a separate column;
there is no gap floor, negative-eigenvalue repair, discarded observation, or
primary raw mean/MAE/RMSE. Exact Wishart sampling uses df=n−1 and is validated
against centered Gaussian observations before pilot execution.

All 32 automated tests passed, with zero errors, failures or skips.
The suite combines the eight existing mathematical tests with 24 E3 checks.
Four Gaussian/Wishart comparisons used 3,000 samples per method, with observed
risk differences 1.19–2.01 combined MCSE (6-MCSE validation limit), alongside
selected theoretical first/second covariance moments. All mandatory projector,
spectrum, normalization, seed, statistical interval, matched bootstrap,
capped-tie, and failure-recording validations passed. A symmetric target-test
fixture was corrected during development because it did not distinguish
population risk from realized loss; this was a test-data correction, not a
change in estimator definition. The later skewed fixture passes.

The smoke test passed all 160 records across four contrasting conditions,
40 each, with exact repeated-seed reproduction and lossless CSV round trip.
Seven raw predictions above 1 were retained. Conservative relative-frequency
intervals were added after smoke review to avoid degenerate bootstrap [0,0]
limits for unobserved events; all tests and smoke were rerun before the pilot.

The full pilot has **84 configurations × 1,000 = 84,000 records**, exactly.
There are zero recorded numerical failures, non-finite observations, or
loss-range violations. Actual losses range from 0.00066829 to 0.58303775.
All 19 postflight checks passed: count/grid integrity, unique IDs,
unchanged protocol, source and artifact hashes, validation/smoke provenance,
independent summary and oracle-classification formulas, exact reconstruction
of two full pilot streams (including the largest extreme), secondary-output
completeness, and identical regenerated numerical tables/report text.

## Primary exploratory findings and uncertainty

- Oracle: 42 within tolerance, 36 outside, 6 uncertain. Plug-in: 12 typically
  within, 68 outside, 4 uncertain. There are 26 oracle-within/plug-in-outside
  configurations. These are pointwise descriptive labels, not confirmed regions.
- Small gaps can fail even at the largest sample size. Fourteen conditions
  have oracle B>1; at n=50, δ=.05, t=0, B=17.13853 while R̂=0.23823
  (95% MC CI 0.23208–0.24439).
- At n=1000, δ=.2, t=0, R̂=0.06218 (0.05728–0.06709), B=0.05821,
  oracle relative error 6.4%, but median relative absolute plug-in error 58.6%
  (matched-bootstrap CI 55.6%–61.6%). Relative underprediction is 41.0%
  (conservative joint-band CI 32.1%–50.3%).
- Among conditions with risk CI above 0.10, serious underprediction peaks at
  44.5% (Wilson CI 41.4%–47.6%) for n=2000, δ=.05, t=0. Here R̂=0.12479
  (0.11792–0.13165), and the raw 95th/5th percentile ratio is 47.07.
- There are 3,032 raw predictions above 1. The maximum is 8,527.07189 at
  n=500, δ=.05, t=.8. Its condition's 99th percentile is 5.07314
  (order-statistic CI 2.40051–14.0815). No finite pilot maximum resolves
  the raw estimator's theoretical infinite mean.

See `pilot_report.md` for the full reliability table and uncertainty methods.
`contrasts.csv` has 42 matched-condition independent-stream comparisons;
`ranking.csv` has 12 seven-gap correlations; `calibration.csv` has 60 groups,
with every fixed-(n,t) stratum accounting for all 7,000 observations.

## Nonmonotone observations required by H1

The estimated oracle relative error increases along 12 adjacent-gap
comparisons in 8 fixed-(n,t) slices. These are descriptive pilot
reversals, not tests of a population trend. The displayed intervals are pointwise
transforms of the existing 95% risk CIs, not confidence intervals for the
between-gap difference; overlapping intervals leave the reversal unresolved.
No hypotheses or thresholds are changed.

| n | t | Adjacent gaps | Earlier relative error [MC interval] | Later relative error [MC interval] |
|---:|---:|---|---:|---:|
| 100 | 0.8 | 0.75 → 1 | 2.23% [0.00%, 5.59%] | 8.02% [4.59%, 11.21%] |
| 250 | 0 | 0.35 → 0.5 | 5.49% [0.00%, 12.01%] | 13.23% [7.36%, 18.40%] |
| 250 | 0.8 | 0.35 → 0.5 | 2.54% [0.00%, 7.43%] | 6.71% [2.13%, 10.88%] |
| 250 | 0.8 | 0.5 → 0.75 | 6.71% [2.13%, 10.88%] | 11.20% [6.89%, 15.12%] |
| 500 | 0.8 | 0.5 → 0.75 | 5.68% [0.64%, 10.24%] | 6.19% [2.62%, 9.52%] |
| 1000 | 0 | 0.2 → 0.35 | 6.39% [0.00%, 13.23%] | 14.46% [6.54%, 21.15%] |
| 1000 | 0 | 0.75 → 1 | 0.01% [0.00%, 4.08%] | 0.17% [0.00%, 3.51%] |
| 1000 | 0.8 | 0.2 → 0.35 | 6.11% [0.00%, 12.12%] | 13.42% [6.98%, 19.03%] |
| 1000 | 0.8 | 0.5 → 0.75 | 3.81% [0.00%, 8.29%] | 4.76% [1.34%, 7.96%] |
| 2000 | 0 | 0.35 → 0.5 | 1.56% [0.00%, 8.26%] | 6.69% [0.49%, 12.16%] |
| 2000 | 0 | 0.75 → 1 | 0.32% [0.00%, 4.08%] | 1.13% [0.00%, 4.52%] |
| 2000 | 0.8 | 0.5 → 0.75 | 3.11% [0.00%, 7.20%] | 4.29% [1.11%, 7.27%] |

## Remaining concerns before confirmation

No unresolved implementation or numerical failure was found. Upper-tail
quantiles at M=1000 remain imprecise. All intervals are pointwise, with no
multiplicity correction. Student and percentile-bootstrap coverage is approximate;
plug-in tolerance margins are nonsmooth statistics. The relative-event CI uses
a conservative DKW/risk-band union rather than relying on a potentially
degenerate percentile bootstrap, while retaining bootstrap limits separately.
Calibration intervals condition on observed bin edges. Seven-configuration
ranking intervals can be narrow or degenerate when bootstrap rankings do not
change; they do not prove a perfect population ranking. Within-sample loss
associations and the additional display of nonmonotone reversals are exploratory.
Review these choices and near-threshold labels before freezing confirmation
analysis. Reserved confirmation namespace 2 has not generated any observations.

## Exact reproduction

To reproduce from a clean checkout, apply the included implementation-only patch
at the recorded base commit, then run:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-e3.txt
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py validate
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py smoke
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py pilot
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py analyze --output results/e3/pilot
```

To regenerate only the supplied results, use the last command in a repository
checkout containing the supplied source and saved results. To recheck the
postflight audit, run:

```sh
OPENBLAS_NUM_THREADS=1 python results/e3/postflight_audit.py
```

Do not overwrite supplied simulation manifests. Repeating the pilot in the same
checkout uses `pilot --output results/e3/pilot-repeat`, after validating and
retaining the default successful smoke. The README describes these safeguards.
Results are tied to the exact recorded software and source, with separate
smoke/pilot/bootstrap streams. No confirmation command is provided.
