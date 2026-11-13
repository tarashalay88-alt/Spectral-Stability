# E3 exploratory pilot report

Run `pilot-e5b01b6e772e-2026111303`: 84 configurations, 84,000 replications; 1000–1000 per configuration. No confirmation was run.

## Mathematical facts and implementation

The target is population risk R=E[L], not the loss of one sample. The Gaussian centered covariance has Wishart law with n−1 degrees of freedom. The pilot uses that exact law after Gaussian/Wishart moment and risk validation. Loss is the squared projector distance divided by 2k; the oracle and raw plug-in use the frozen formula divided by k(n−1). Raw values have no cap or gap floor; the separate capped estimate is min(1,max(0,B̂)). No eigenvalues are repaired. The spectrum, trace 23, leading eigenvalue 6, and fixed top-3 population subspace are tested across the full frozen grid. The old clipped loss helper is used only as an independent validation reference; E3 retains the unclipped computed loss.

The protocol states that the raw plug-in has an infinite mathematical mean. A finite empirical sample neither verifies nor contradicts this fact. Medians, quantiles, and event frequencies are primary; capped MAE/RMSE use R̂.

## Pilot observations

Oracle labels: {'within': 42, 'outside': 36, 'uncertain': 6}. Typical plug-in labels: {'outside': 68, 'within': 12, 'uncertain': 4}. There are 26 configurations with oracle within tolerance and plug-in typically outside tolerance. Labels use max(0.002,0.20R) and include uncertainty.

Oracle relative error ranges from 0.0% (n=1000, δ=0.75, t=0) to 7094.1% (n=50, δ=0.05, t=0). 14 configurations have B>1. These bounds on R are a mathematical breakdown diagnostic; the observed error maps provide the configuration-specific comparisons.

| Configuration | R̂ [95% MC CI] | B | Oracle relative error | Median relative absolute plug-in error | Oracle / plug-in label |
|---|---:|---:|---:|---:|---|
| c000: n=50, δ=0.05, t=0 | 0.23823 [0.23208, 0.24439] | 17.13853 | 7094.1% | 40.9% | outside / outside |
| c066: n=1000, δ=0.75, t=0 | 0.00872 [0.00838, 0.00906] | 0.00872 | 0.0% | 13.4% | within / within |
| c082: n=2000, δ=1, t=0 | 0.00333 [0.00323, 0.00344] | 0.00337 | 1.1% | 6.3% | within / within |
| c043: n=500, δ=0.05, t=0.8 | 0.16437 [0.15761, 0.17112] | 1.70269 | 935.9% | 57.6% | outside / outside |
| c010: n=50, δ=0.75, t=0 | 0.16311 [0.15777, 0.16846] | 0.17786 | 9.0% | 36.4% | within / outside |
| c070: n=2000, δ=0.05, t=0 | 0.12479 [0.11792, 0.13165] | 0.42010 | 236.7% | 70.2% | outside / outside |

| n | t | Gaps with oracle within tolerance | Oracle uncertain count | Plug-in typically within / total |
|---:|---:|---|---:|---:|
| 50 | 0 | 0.75, 1 | 0 | 0 / 7 |
| 50 | 0.8 | 1 | 1 | 0 / 7 |
| 100 | 0 | 0.5, 0.75, 1 | 0 | 0 / 7 |
| 100 | 0.8 | 0.75, 1 | 1 | 0 / 7 |
| 250 | 0 | 0.35, 0.5, 0.75, 1 | 0 | 0 / 7 |
| 250 | 0.8 | 0.35, 0.5, 0.75, 1 | 0 | 0 / 7 |
| 500 | 0 | 0.35, 0.5, 0.75, 1 | 1 | 1 / 7 |
| 500 | 0.8 | 0.35, 0.5, 0.75, 1 | 0 | 1 / 7 |
| 1000 | 0 | 0.2, 0.5, 0.75, 1 | 1 | 2 / 7 |
| 1000 | 0.8 | 0.2, 0.35, 0.5, 0.75, 1 | 0 | 2 / 7 |
| 2000 | 0 | 0.2, 0.35, 0.5, 0.75, 1 | 1 | 3 / 7 |
| 2000 | 0.8 | 0.35, 0.5, 0.75, 1 | 1 | 3 / 7 |

The largest raw 95th/5th percentile ratio is 47.07 at n=2000, δ=0.05, t=0 (5th=0.015655, 95th=0.73683). This measures within-condition prediction spread, separately from typical error.

The largest retained raw prediction is 8527.07 at n=500, δ=0.05, t=0.8; the raw 99th percentile there is 5.07314 [95% CI 2.40051, 14.0815]. Across all observations, 3,032 raw values exceed 1. Prediction-spread plots display all six frozen quantiles (including the median).

Even among oracle-within-tolerance configurations, relative underprediction reaches 41.0% [95% CI 32.1%, 50.3%] at n=1000, δ=0.2, t=0. Its median relative absolute plug-in error is 58.6% [95% CI 55.6%, 61.6%].

Among 42 configurations whose risk CI is entirely above 0.10, the largest serious-underprediction frequency is 44.5% [95% Wilson CI 41.4%, 47.6%] at n=2000, δ=0.05, t=0. Relative underprediction there is 54.0% [conservative joint-band CI 46.1%, 60.5%].

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
