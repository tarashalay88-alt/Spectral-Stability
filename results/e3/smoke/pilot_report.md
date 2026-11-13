# E3 exploratory smoke report

Run `smoke-e5b01b6e772e-2026111303`: 4 configurations, 160 replications; 40–40 per configuration. No confirmation was run.

## Mathematical facts and implementation

The target is population risk R=E[L], not the loss of one sample. The Gaussian centered covariance has Wishart law with n−1 degrees of freedom. The pilot uses that exact law after Gaussian/Wishart moment and risk validation. Loss is the squared projector distance divided by 2k; the oracle and raw plug-in use the frozen formula divided by k(n−1). Raw values have no cap or gap floor; the separate capped estimate is min(1,max(0,B̂)). No eigenvalues are repaired. The spectrum, trace 23, leading eigenvalue 6, and fixed top-3 population subspace are tested across the full frozen grid. The old clipped loss helper is used only as an independent validation reference; E3 retains the unclipped computed loss.

The protocol states that the raw plug-in has an infinite mathematical mean. A finite empirical sample neither verifies nor contradicts this fact. Medians, quantiles, and event frequencies are primary; capped MAE/RMSE use R̂.

## Pilot observations

Oracle labels: {'outside': 2, 'within': 2}. Typical plug-in labels: {'outside': 2, 'within': 2}. There are 0 configurations with oracle within tolerance and plug-in typically outside tolerance. Labels use max(0.002,0.20R) and include uncertainty.

Oracle relative error ranges from 0.3% (n=2000, δ=1, t=0.8) to 6866.6% (n=50, δ=0.05, t=0). 2 configurations have B>1. These bounds on R are a mathematical breakdown diagnostic; the observed error maps provide the configuration-specific comparisons.

| Configuration | R̂ [95% MC CI] | B | Oracle relative error | Median relative absolute plug-in error | Oracle / plug-in label |
|---|---:|---:|---:|---:|---|
| c000: n=50, δ=0.05, t=0 | 0.24601 [0.21228, 0.27973] | 17.13853 | 6866.6% | 41.1% | outside / outside |
| c083: n=2000, δ=1, t=0.8 | 0.00546 [0.00484, 0.00608] | 0.00547 | 0.3% | 10.8% | within / within |
| c082: n=2000, δ=1, t=0 | 0.00361 [0.00300, 0.00421] | 0.00337 | 6.5% | 7.4% | within / within |
| c001: n=50, δ=0.05, t=0.8 | 0.27149 [0.23852, 0.30445] | 17.33964 | 6286.9% | 44.2% | outside / outside |

| n | t | Gaps with oracle within tolerance | Oracle uncertain count | Plug-in typically within / total |
|---:|---:|---|---:|---:|
| 50 | 0 | none | 0 | 0 / 1 |
| 50 | 0.8 | none | 0 | 0 / 1 |
| 2000 | 0 | 1 | 0 | 1 / 1 |
| 2000 | 0.8 | 1 | 0 | 1 / 1 |

The largest raw 95th/5th percentile ratio is 17.47 at n=50, δ=0.05, t=0.8 (5th=0.11057, 95th=1.9318). This measures within-condition prediction spread, separately from typical error.

The largest retained raw prediction is 4.69717 at n=50, δ=0.05, t=0.8; the raw 99th percentile there is 3.95434 [95% CI 1.88653, inf]. Across all observations, 7 raw values exceed 1. Prediction-spread plots display all six frozen quantiles (including the median).

Even among oracle-within-tolerance configurations, relative underprediction reaches 0.0% [95% CI 0.0%, 23.4%] at n=2000, δ=1, t=0. Its median relative absolute plug-in error is 7.4% [95% CI 4.5%, 18.4%].

Among 2 configurations whose risk CI is entirely above 0.10, the largest serious-underprediction frequency is 0.0% [95% Wilson CI 0.0%, 8.8%] at n=50, δ=0.05, t=0. Relative underprediction there is 2.5% [conservative joint-band CI 0.0%, 45.9%].

## Uncertainty and secondary analyses

All intervals are pointwise 95%; there is no multiplicity correction or claim of a universal transition boundary. R̂ uses the sample SD/√M and a Student-t interval. Raw quantiles use conservative binomial order-statistic intervals; frequencies with fixed thresholds use Wilson intervals. An infinite quantile CI endpoint in a small smoke sample means insufficient distribution-free tail information, not a non-finite simulated observation. Joint percentile bootstrap intervals use 2000 resamples, resample complete (L,B̂,B̂cap,gap) records, and recompute R̂ in every draw. Plug-in labels use the bootstrap CI for median(|B̂−R̂|)−max(0.002,0.20R̂). Oracle labels check the whole risk interval, including the tolerance breakpoint.

Relative-underprediction intervals use a conservative union of a 97.5% Student risk interval and a 97.5% uniform DKW band for the raw-prediction CDF. This accounts for the estimated threshold and remains nondegenerate when no event is observed, without assuming independence of loss and prediction. Coverage inherits the Student interval approximation; the CDF component is distribution-free. Matched-bootstrap intervals are retained separately as `relative_underprediction_bootstrap_low/high` for comparison.

Streams are independent across configurations, so matched (n,δ) background contrasts resample the two configurations independently. `contrasts.csv` reports risk, median-prediction, and median-absolute-error differences. `ranking.csv` compares only seven gaps at fixed (n,t), with bootstrap intervals. `calibration.csv` pools those seven gaps with equal simulated weights, merges tied quantile edges, and uses a stratified complete-record bootstrap that recomputes each reference risk. Its intervals are conditional on the observed displayed bin boundaries. Secondary tables are omitted from smoke runs with incomplete seven-gap strata. Within-sample loss/prediction and gap/prediction Spearman associations in `summary.csv` are explicitly exploratory, without confidence or causal claims.

## Exploratory interpretation and concerns before confirmation

Inspect fixed-(n,t) slices of the continuous maps: pooled results cannot isolate a gap effect. Small typical error does not imply narrow prediction spread. The robust decomposition columns separate median(B̂−B) from B−R̂; cancellation between them is not evidence for oracle accuracy. No strict monotonicity or new theorem is asserted. Ambiguous labels remain uncertain.

Numerical failures: 0; non-finite replication records: 0; loss-range violations beyond 1e−12: 0. Extreme finite raw observations are preserved, not discarded. `summary.csv` includes minima of the observed gaps. Zero failures do not remove heavy-tail risk.

Before confirmation, review Monte Carlo precision near the frozen labels, the stability of upper-tail quantiles at M=1000, the approximate coverage of percentile bootstrap intervals for nonsmooth median-error labels, and calibration intervals conditional on estimated bin boundaries. Treat seven-point rank correlations cautiously. Eligibility of serious underprediction is uncertain when the risk CI crosses 0.10; frequencies are retained for all configurations but the serious-underprediction map shows only established eligibility. Review this analysis implementation and its software choices before freezing confirmation analysis. No amendment to hypotheses, grid, definitions, or thresholds is proposed or applied. Confirmation namespace 2 is reserved and unused.

## Reproduction and provenance

Base repository commit: `afa3dcbdb54909d82193c4a9d6f78d785d5b62db`. The implementation is identified exactly by per-file SHA-256 hashes and a source fingerprint in `run_manifest.json`, which also records the dirty Git state, all seed keys, parameters, Python/library versions, environment, numerical counts, and artifact hashes. `validation.json` and `smoke/run_manifest.json` are linked by hash from the pilot manifest. `figure_index.json` links each figure to this run.

Run from the repository root (see `experiments/e3/README.md`):

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-e3.txt
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py validate
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py smoke
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py pilot
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py analyze --output results/e3/pilot
```

Simulation commands refuse to overwrite existing run manifests; use a fresh `--output` directory for a repeat. `analyze` regenerates tables, figures, and this report from the retained replication CSV and its recorded configuration. No confirmatory simulation entry point is implemented.
