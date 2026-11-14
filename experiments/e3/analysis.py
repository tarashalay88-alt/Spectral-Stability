"""E3-A/E3-B and prespecified secondary outputs from retained matched rows."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from PCA_stability import first_order_benchmark
from spectral_stability.e3_statistics import (
    METRICS, QUANTILES, bootstrap_metrics, calibration_edges, calibration_groups,
    joint_metrics, mean_interval, oracle_label, percentile_interval,
    quantile_interval, relative_frequency_interval, tolerance, wilson_interval)
from .design import eigenvalues, seed_sequence
from .simulation import VALUE_COLUMNS


def summarize(records, config, stage):
    if records.failure.any() or not np.isfinite(records[VALUE_COLUMNS]).all().all():
        raise ValueError('Numerical failures retained; investigate before inference')
    summaries, bootstrap = [], {}
    for config_id, group in records.groupby('config_id', sort=True):
        row = {key:group.iloc[0][key] for key in
               ('config_id', 'config_index', 'n', 'delta', 't')}
        row['M'] = len(group)
        values = group[VALUE_COLUMNS].to_numpy()
        raw, loss = values[:, 1], values[:, 0]
        benchmark = first_order_benchmark(np.diag(eigenvalues(row['delta'], row['t'])),
                                          int(row['n']), config['k'])
        risk, se, lower, upper = mean_interval(loss)
        row.update(risk=risk, risk_mcse=se, risk_low=lower, risk_high=upper,
                   oracle=benchmark, oracle_signed_error=benchmark-risk,
                   oracle_absolute_error=abs(benchmark-risk),
                   oracle_relative_error=abs(benchmark-risk)/risk,
                   tolerance=float(tolerance(risk)),
                   oracle_label=oracle_label(benchmark, lower, upper))
        estimates = joint_metrics(values, benchmark)
        rng = np.random.default_rng(seed_sequence(config, stage, int(row['config_index']), 1))
        draws = bootstrap_metrics(values, benchmark, rng, config['bootstrap_repetitions'])
        bootstrap[config_id] = draws
        lo, hi = percentile_interval(draws)
        for i, name in enumerate(METRICS):
            if name == 'risk':  # primary risk interval remains Student-t
                continue
            row[name], row[name+'_low'], row[name+'_high'] = estimates[i], lo[i], hi[i]
        for endpoint in ('low','high'):
            row['relative_underprediction_bootstrap_'+endpoint] = row['relative_underprediction_'+endpoint]
        (row['relative_underprediction_low'], row['relative_underprediction_high']) = relative_frequency_interval(loss,raw)
        row['plugin_label'] = ('within' if row['plugin_tolerance_margin_high'] <= 0 else
                               'outside' if row['plugin_tolerance_margin_low'] > 0 else
                               'uncertain')
        for probability in QUANTILES:
            name = f'raw_q{int(100*probability):02d}'
            row[name] = np.quantile(raw, probability)
            row[name+'_low'], row[name+'_high'] = quantile_interval(raw, probability)
        for name, events in [('raw_exceeds_one', raw > 1), ('serious_underprediction', raw <= .05)]:
            row[name] = events.mean()
            row[name+'_count'] = int(events.sum())
            row[name+'_low'], row[name+'_high'] = wilson_interval(events.sum(), len(raw))
        row['serious_eligibility'] = ('established' if lower > .10 else
                                      'below_threshold' if upper < .10 else 'uncertain')
        row['maximum_raw'] = raw.max()
        row['minimum_sample_gap'] = values[:, 3].min()
        row['failure_count'] = int(group.failure.sum())
        row['loss_raw_spearman'] = spearmanr(loss, raw).statistic
        row['gap_raw_spearman'] = spearmanr(values[:, 3], raw).statistic
        # Robust decomposition: median(raw-B)+(B-R)=median(raw)-R.
        row['decomposition_residual'] = (row['median_sample_eigenvalue_error']
                                       +row['oracle_signed_error']-row['median_signed_error'])
        summaries.append(row)
    return pd.DataFrame(summaries), bootstrap


def contrasts_and_ranking(summary, bootstrap):
    contrasts, ranking = [], []
    for (n, delta), group in summary.groupby(['n', 'delta']):
        if len(group) != 2:
            continue
        base, alternative = group.sort_values('t').to_dict('records')
        row = dict(n=n, delta=delta, reference_t=0., alternative_t=.8,
                   design='independent configuration streams')
        for name in ('risk', 'median_prediction', 'median_absolute_error'):
            i = METRICS.index(name)
            row[name+'_difference'] = alternative[name]-base[name]
            lo, hi = percentile_interval(bootstrap[alternative['config_id']][:, i]
                                         -bootstrap[base['config_id']][:, i])
            row[name+'_difference_low'], row[name+'_difference_high'] = lo, hi
        row['oracle_difference'] = alternative['oracle']-base['oracle']
        contrasts.append(row)
    for (n, t), group in summary.groupby(['n', 't']):
        if len(group) != 7:
            continue
        group = group.sort_values('delta')
        draws = np.stack([bootstrap[x] for x in group.config_id])
        a = rankdata(draws[:, :, METRICS.index('risk')], axis=0)
        b = rankdata(draws[:, :, METRICS.index('median_prediction')], axis=0)
        a, b = a-a.mean(axis=0), b-b.mean(axis=0)
        correlations = np.sum(a*b, axis=0)/np.sqrt(np.sum(a*a, axis=0)*np.sum(b*b, axis=0))
        lo, hi = percentile_interval(correlations)
        ranking.append(dict(n=n, t=t, configurations=7,
                            spearman=spearmanr(group.risk, group.median_prediction).statistic,
                            ci_low=lo, ci_high=hi))
    return pd.DataFrame(contrasts), pd.DataFrame(ranking)


def calibration(records, summary, config, stage):
    """Fixed observed quantile bins; stratified matched-row bootstrap.

    Intervals are conditional on the displayed bin boundaries, and include
    uncertainty in every reference R. No tie is split, even in resampling.
    """
    result = []
    references = summary.set_index('config_id').risk
    for (n, t), frame in records.groupby(['n', 't']):
        if frame.config_id.nunique() != 7:
            continue
        edges = calibration_edges(frame.plugin_capped.to_numpy())
        labels = calibration_groups(frame.plugin_capped.to_numpy(), edges)
        groups = len(edges)-1
        reference = frame.config_id.map(references).to_numpy()
        observed = []
        for label in range(groups):
            mask = labels == label
            observed.append(dict(n=n, t=t, group=label, groups_after_ties=groups,
                                 lower_edge=edges[label], upper_edge=edges[label+1],
                                 count=int(mask.sum()), mean_prediction=frame.plugin_capped[mask].mean(),
                                 mean_reference_risk=reference[mask].mean(),
                                 capped_ties=int((frame.plugin_capped[mask] == 1).sum())))
        first_index = int(frame.config_index.min())
        rng = np.random.default_rng(seed_sequence(config, stage, first_index, 2))
        draws = []
        strata = [g[VALUE_COLUMNS].to_numpy() for _,g in frame.groupby('config_id', sort=True)]
        for start in range(0, config['bootstrap_repetitions'], 100):
            size = min(100, config['bootstrap_repetitions']-start)
            predictions, risks = [], []
            for values in strata:
                index = rng.integers(0, len(values), (size, len(values)))
                predictions.append(values[index, 2])
                risk = values[index, 0].mean(axis=1)
                risks.append(np.broadcast_to(risk[:, None], index.shape))
            predictions, risks = np.concatenate(predictions, axis=1), np.concatenate(risks, axis=1)
            group_ids = calibration_groups(predictions, edges)
            statistics = []
            for label in range(groups):
                mask = group_ids == label
                counts = mask.sum(axis=1)
                if np.any(counts == 0):
                    raise ValueError('Empty calibration bootstrap group; investigate')
                p = (predictions*mask).sum(axis=1)/counts
                r = (risks*mask).sum(axis=1)/counts
                statistics.append(np.stack([p,r,p-r], axis=1))
            draws.append(np.stack(statistics, axis=1))
        lo, hi = percentile_interval(np.concatenate(draws))
        for label, row in enumerate(observed):
            row['prediction_minus_reference'] = row['mean_prediction']-row['mean_reference_risk']
            for i,name in enumerate(('mean_prediction','mean_reference_risk','prediction_minus_reference')):
                row[name+'_low'], row[name+'_high'] = lo[label,i], hi[label,i]
            result.append(row)
    return pd.DataFrame(result)
