# E3 pilot

`E3.md` is the frozen protocol and is unchanged. This implementation runs the
84-condition, 1,000-replication exploratory pilot only. There is no confirmation
command. Python 3.12 was used; exact libraries are in `requirements-e3.txt`.

From the repository root:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-e3.txt
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py validate
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py smoke
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py pilot
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py analyze --output results/e3/pilot
```

`validate` runs all legacy mathematical tests and the E3 suite, records its log,
source hashes and pass status in `results/e3/validation.*`. It includes projector
properties, independent angle formulas, hand calculations, the full frozen
spectrum, Gaussian/Wishart moments and risk comparisons, interval coverage,
bootstrap dependence, capped ties, seed separation, and failure recording.
`smoke` requires current successful validation and runs 40 replications for each
of (n=50,delta=.05) and (n=2000,delta=1), for both backgrounds. It repeats the
draws to verify exact seed reproduction, verifies the CSV round trip, and runs
the full analysis. `pilot` requires successful current validation and default
smoke artifacts. Changing source, tests, configuration or protocol makes those
gates stale. Any recorded numerical failure blocks statistical analysis.

Simulation commands refuse existing run manifests. To repeat without overwriting:

```sh
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py smoke --output results/e3/smoke-repeat
OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py pilot --output results/e3/pilot-repeat
```

A repeated pilot still checks the default smoke gate. `analyze` reads the saved
CSV, checks its hash and the source hashes, then regenerates all analysis outputs
without further simulation. Use `--output` to choose the saved run. Gate artifacts
are in `results/e3/`; run-specific files are in its `smoke/` and `pilot/` directories.

## Mathematical and statistical choices

The reusable core is `src/spectral_stability/e3_math.py`; inference is in
`e3_statistics.py`. Existing Gaussian observations, centered covariance and
oracle benchmark utilities are reused. A single symmetric eigendecomposition
provides the full descending sample spectrum and principal subspace. The
Frobenius projector loss is retained without clipping, with range tolerance
1e-12. The older clipped loss helper is an independent test reference only.
The raw benchmark has neither gap floor nor cap; the explicitly separate capped
column alone is bounded. Negative sample eigenvalues are diagnosed, never repaired.

SciPy Bartlett Wishart sampling, divided by n-1, has exactly the centered Gaussian
sample covariance law. Four contrasting configurations are compared against
direct centered Gaussian sampling in validation, using 3,000 replications per
method and 6-MCSE tolerances. Selected diagonal/off-diagonal first and second
moments use their theoretical Wishart values. Rotated covariances and small
singular samples are tested separately. Reusable sampling allows general d,n;
df<d uses direct observations, and general d,k,n measurements are tested.

PCG64 uses `SeedSequence(2026111303, spawn_key=(stage,config_index,purpose))`.
Stages 0/1/2 are smoke/pilot/reserved confirmation; purposes 0/1/2 are simulation,
joint bootstrap and secondary calibration. Configurations use the full grid
order n, delta, t; smoke retains those grid IDs. No common random numbers are
used across configurations. Secondary bootstrap contrasts and ranking combine
independently resampled configurations. Every stream is deterministically
recoverable from the configuration, stage and purpose; confirmation is unused.

R uses s_L/sqrt(M) with a Student-t 95% interval. Oracle classification evaluates
the tolerance over the entire risk CI, including piecewise-linear breakpoints.
All raw quantiles have binomial order-statistic CIs (conservative with ties).
These can have infinite endpoints in small smoke samples. Fixed-threshold
frequencies have Wilson CIs. Other plug-in measures and labels use 2,000
complete-record percentile bootstrap draws, recomputing risk in every draw.
All intervals are pointwise, approximate where stated, and not adjusted for
84 comparisons. No primary raw sample mean, raw MAE or raw RMSE is calculated.

The primary relative-underprediction CI combines a 97.5% Student risk interval
with a 97.5% uniform DKW band for the raw-prediction CDF (union bound). It accounts
for the estimated risk threshold and is conservative without an independence
assumption; risk coverage remains approximate. Complete-record percentile
bootstrap limits are retained separately with `_bootstrap_low/high` suffixes.
The joint-band interval avoids the degenerate [0,0] bootstrap interval possible
with zero observed events. These choices were implemented before the pilot.

Summary columns `median_absolute_error` and `relative_median_absolute_error`
mean median(|raw-Rhat|) and its ratio to Rhat, respectively. They do not refer
to the median absolute deviation about the sample median or to realized loss.
`median_sample_eigenvalue_error` plus `oracle_signed_error` equals
`median_signed_error`; the residual is checked. The two tolerance labels are
kept separate, including uncertainty. Serious-underprediction frequencies are
saved for all conditions, with risk-based eligibility recorded separately;
the primary serious map displays only risk CIs entirely above 0.10.

Secondary ranking uses seven medians/risks per fixed (n,t); bootstrap CIs reflect
independent configuration streams. Calibration pools the seven gaps equally,
merges tied quantile cutpoints, and never splits equal capped values. Its
stratified matched bootstrap includes reference-risk uncertainty and is
conditional on the observed bin boundaries. Within-sample Spearman associations
of raw predictions with loss and gap are exploratory descriptive columns only.

## Outputs and audit

`replications.csv` stores every matched value, condition, replication ID, stream
key and failure flag. `summary.csv` includes all primary point estimates and CIs;
`contrasts.csv`, `ranking.csv`, `calibration.csv` contain secondary analyses.
`figures/` holds PNG and SVG plots with run and manifest links. Missing smoke
cells stay blank; plots do not extrapolate. `pilot_report.md` separates facts,
observations, interpretation and remaining uncertainties. The same filename is
used for the smoke report, with its stage clearly stated.

`run_manifest.json` records the base Git commit, working-tree state, exact source
hashes, protocol hash, source fingerprint, complete configuration and seed keys,
Python/library/BLAS/environment details, replication counts, simulation runtime,
numerical failures, validation/smoke hashes and analysis-artifact hashes.
`figure_index.json` links figures to that manifest and source fingerprint.
Results depend on the recorded NumPy/SciPy RNG implementation; reproducibility
of numerical bits across different BLAS libraries/platforms is not asserted.

Validation log and manifest are the authoritative gate records. The repository
top-level `test.py` is a single-draw demonstration, not an automated test suite;
it is intentionally not treated as a Monte Carlo validation. No changes were
made to E1/E2 or to the frozen E3 protocol.
