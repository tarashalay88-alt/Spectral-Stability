"""Concise exploratory report, regenerated from data rather than hand edits."""
from pathlib import Path


def condition(row):
    return f"n={row.n}, δ={row.delta:g}, t={row.t:g}"


def make_report(summary, records, manifest, output):
    oracle = summary.oracle_label.value_counts().to_dict()
    plugin = summary.plugin_label.value_counts().to_dict()
    mismatch = summary[(summary.oracle_label=='within') & (summary.plugin_label=='outside')]
    eligible = summary[summary.serious_eligibility=='established']
    peak = summary.loc[summary.maximum_raw.idxmax()]
    worst = summary.loc[summary.oracle_relative_error.idxmax()]
    best = summary.loc[summary.oracle_relative_error.idxmin()]
    typical_best = summary.loc[summary.relative_median_absolute_error.idxmin()]
    target = eligible.loc[eligible.serious_underprediction.idxmax()] if len(eligible) else None
    reliable = summary[summary.oracle_label=='within']
    reliability_rows = []
    for (n,t), group in summary.groupby(['n','t']):
        gaps = ', '.join(f'{x:g}' for x in group.loc[group.oracle_label=='within','delta']) or 'none'
        reliability_rows.append(f'| {n} | {t:g} | {gaps} | '
            f"{int((group.oracle_label=='uncertain').sum())} | "
            f"{int((group.plugin_label=='within').sum())} / {len(group)} |")
    widest = summary.loc[(summary.raw_q95/summary.raw_q05).idxmax()]
    examples = [worst, best, typical_best, peak]
    if len(mismatch):
        examples.append(mismatch.iloc[0])
    if target is not None:
        examples.append(target)
    seen, rows = set(), []
    for row in examples:
        if row.config_id in seen:
            continue
        seen.add(row.config_id)
        rows.append(f"| {row.config_id}: {condition(row)} | {row.risk:.5f} "
                    f"[{row.risk_low:.5f}, {row.risk_high:.5f}] | {row.oracle:.5f} | "
                    f"{row.oracle_relative_error:.1%} | {row.relative_median_absolute_error:.1%} | "
                    f"{row.oracle_label} / {row.plugin_label} |")
    lines = [f"# E3 exploratory {manifest['stage']} report", '',
             f"Run `{manifest['run_id']}`: {len(summary)} configurations, {len(records):,} replications; "
             f"{summary.M.min()}–{summary.M.max()} per configuration. No confirmation was run.", '',
             '## Mathematical facts and implementation', '',
             'The target is population risk R=E[L], not the loss of one sample. '
             'The Gaussian centered covariance has Wishart law with n−1 degrees of freedom. '
             'The pilot uses that exact law after Gaussian/Wishart moment and risk validation. '
             'Loss is the squared projector distance divided by 2k; the oracle and raw plug-in '
             'use the frozen formula divided by k(n−1). Raw values have no cap or gap floor; '
             'the separate capped estimate is min(1,max(0,B̂)). No eigenvalues are repaired. '
             'The spectrum, trace 23, leading eigenvalue 6, and fixed top-3 population subspace '
             'are tested across the full frozen grid. The old clipped loss helper is used only '
             'as an independent validation reference; E3 retains the unclipped computed loss.', '',
             'The protocol states that the raw plug-in has an infinite mathematical mean. '
             'A finite empirical sample neither verifies nor contradicts this fact. '
             'Medians, quantiles, and event frequencies are primary; capped MAE/RMSE use R̂.', '',
             '## Pilot observations', '',
             f"Oracle labels: {oracle}. Typical plug-in labels: {plugin}. "
             f"There are {len(mismatch)} configurations with oracle within tolerance and plug-in "
             'typically outside tolerance. Labels use max(0.002,0.20R) and include uncertainty.', '',
             f"Oracle relative error ranges from {best.oracle_relative_error:.1%} ({condition(best)}) "
             f"to {worst.oracle_relative_error:.1%} ({condition(worst)}). "
             f"{int((summary.oracle>1).sum())} configurations have B>1. "
             'These bounds on R are a mathematical breakdown diagnostic; the observed error '
             'maps provide the configuration-specific comparisons.', '',
             '| Configuration | R̂ [95% MC CI] | B | Oracle relative error | Median relative absolute plug-in error | Oracle / plug-in label |',
             '|---|---:|---:|---:|---:|---|', *rows, '',
             '| n | t | Gaps with oracle within tolerance | Oracle uncertain count | Plug-in typically within / total |',
             '|---:|---:|---|---:|---:|', *reliability_rows, '',
             f"The largest raw 95th/5th percentile ratio is "
             f"{widest.raw_q95/widest.raw_q05:.2f} at {condition(widest)} "
             f"(5th={widest.raw_q05:.5g}, 95th={widest.raw_q95:.5g}). "
             'This measures within-condition prediction spread, separately from typical error.', '',
             f"The largest retained raw prediction is {peak.maximum_raw:.6g} at {condition(peak)}; "
             f"the raw 99th percentile there is {peak.raw_q99:.6g} "
             f"[95% CI {peak.raw_q99_low:.6g}, {peak.raw_q99_high:.6g}]. "
             f"Across all observations, {int((records.plugin_raw>1).sum()):,} raw values exceed 1. "
             'Prediction-spread plots display all six frozen quantiles (including the median).', '']
    if len(reliable):
        unreliable_sample = reliable.loc[reliable.relative_underprediction.idxmax()]
        lines += [f"Even among oracle-within-tolerance configurations, relative underprediction "
                  f"reaches {unreliable_sample.relative_underprediction:.1%} "
                  f"[95% CI {unreliable_sample.relative_underprediction_low:.1%}, "
                  f"{unreliable_sample.relative_underprediction_high:.1%}] at "
                  f"{condition(unreliable_sample)}. Its median relative absolute plug-in error is "
                  f"{unreliable_sample.relative_median_absolute_error:.1%} "
                  f"[95% CI {unreliable_sample.relative_median_absolute_error_low:.1%}, "
                  f"{unreliable_sample.relative_median_absolute_error_high:.1%}].", '']
    if target is not None:
        lines += [f"Among {len(eligible)} configurations whose risk CI is entirely above 0.10, "
                  f"the largest serious-underprediction frequency is {target.serious_underprediction:.1%} "
                  f"[95% Wilson CI {target.serious_underprediction_low:.1%}, "
                  f"{target.serious_underprediction_high:.1%}] at {condition(target)}. "
                  f"Relative underprediction there is {target.relative_underprediction:.1%} "
                  f"[conservative joint-band CI {target.relative_underprediction_low:.1%}, "
                  f"{target.relative_underprediction_high:.1%}].", '']
    lines += ['## Uncertainty and secondary analyses', '',
              'All intervals are pointwise 95%; there is no multiplicity correction or claim '
              'of a universal transition boundary. R̂ uses the sample SD/√M and a Student-t '
              'interval. Raw quantiles use conservative binomial order-statistic intervals; '
              'frequencies with fixed thresholds use Wilson intervals. An infinite quantile '
              'CI endpoint in a small smoke sample means insufficient distribution-free '
              'tail information, not a non-finite simulated observation. Joint percentile '
              f"bootstrap intervals use {manifest['config']['bootstrap_repetitions']} resamples, "
              'resample complete (L,B̂,B̂cap,gap) records, and recompute R̂ in every draw. '
              'Plug-in labels use the bootstrap CI for median(|B̂−R̂|)−max(0.002,0.20R̂). '
              'Oracle labels check the whole risk interval, including the tolerance breakpoint.', '',
              'Relative-underprediction intervals use a conservative union of a 97.5% '
              'Student risk interval and a 97.5% uniform DKW band for the raw-prediction CDF. '
              'This accounts for the estimated threshold and remains nondegenerate when no '
              'event is observed, without assuming independence of loss and prediction. '
              'Coverage inherits the Student interval approximation; the CDF component is '
              'distribution-free. Matched-bootstrap intervals are retained separately as '
              '`relative_underprediction_bootstrap_low/high` for comparison.', '',
              'Streams are independent across configurations, so matched (n,δ) background '
              'contrasts resample the two configurations independently. `contrasts.csv` '
              'reports risk, median-prediction, and median-absolute-error differences. '
              '`ranking.csv` compares only seven gaps at fixed (n,t), with bootstrap intervals. '
              '`calibration.csv` pools those seven gaps with equal simulated weights, merges '
              'tied quantile edges, and uses a stratified complete-record bootstrap that '
              'recomputes each reference risk. Its intervals are conditional on the observed '
              'displayed bin boundaries. Secondary tables are omitted from smoke runs with '
              'incomplete seven-gap strata. Within-sample loss/prediction and gap/prediction '
              'Spearman associations in `summary.csv` are explicitly exploratory, without '
              'confidence or causal claims.', '',
              '## Exploratory interpretation and concerns before confirmation', '',
              'Inspect fixed-(n,t) slices of the continuous maps: pooled results cannot isolate '
              'a gap effect. Small typical error does not imply narrow prediction spread. '
              'The robust decomposition columns separate median(B̂−B) from B−R̂; cancellation '
              'between them is not evidence for oracle accuracy. No strict monotonicity or '
              'new theorem is asserted. Ambiguous labels remain uncertain.', '',
              f"Numerical failures: {int(records.failure.sum())}; non-finite replication records: "
              f"{int(records.nonfinite.sum())}; loss-range violations beyond 1e−12: "
              f"{int(records.loss_out_of_bounds.sum())}. Extreme finite raw observations are "
              'preserved, not discarded. `summary.csv` includes minima of the observed gaps. '
              'Zero failures do not remove heavy-tail risk.', '',
              'Before confirmation, review Monte Carlo precision near the frozen labels, '
              'the stability of upper-tail quantiles at M=1000, the approximate coverage of '
              'percentile bootstrap intervals for nonsmooth median-error labels, and '
              'calibration intervals conditional on estimated bin boundaries. Treat seven-point '
              'rank correlations cautiously. Eligibility of serious underprediction is uncertain '
              'when the risk CI crosses 0.10; frequencies are retained for all configurations '
              'but the serious-underprediction map shows only established eligibility. '
              'Review this analysis implementation and its software choices before freezing '
              'confirmation analysis. No amendment to hypotheses, grid, definitions, or '
              'thresholds is proposed or applied. Confirmation namespace 2 is reserved and unused.', '',
              '## Reproduction and provenance', '',
              f"Base repository commit: `{manifest['git_base_commit']}`. The implementation is "
              'identified exactly by per-file SHA-256 hashes and a source fingerprint in '
              '`run_manifest.json`, which also records the dirty Git state, all seed keys, '
              'parameters, Python/library versions, environment, numerical counts, and '
              'artifact hashes. `validation.json` and `smoke/run_manifest.json` are linked '
              'by hash from the pilot manifest. `figure_index.json` links each figure to this run.', '',
              'Run from the repository root (see `experiments/e3/README.md`):', '',
              '```sh', 'python -m venv .venv', '. .venv/bin/activate',
              'python -m pip install -r requirements-e3.txt',
              'OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py validate',
              'OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py smoke',
              'OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py pilot',
              'OPENBLAS_NUM_THREADS=1 python experiments/e3/run.py analyze --output results/e3/pilot',
              '```', '',
              'Simulation commands refuse to overwrite existing run manifests; use a fresh '
              '`--output` directory for a repeat. `analyze` regenerates tables, figures, and '
              'this report from the retained replication CSV and its recorded configuration. '
              'No confirmatory simulation entry point is implemented.', '']
    Path(output).write_text('\n'.join(lines))
