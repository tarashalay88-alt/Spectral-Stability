"""Independent-replication Monte Carlo inference; no raw mean/MAE/RMSE.

All intervals are pointwise 95%, not simultaneous confidence statements.
The joint nonparametric bootstrap uses one index matrix for complete records.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import binom, norm, t

QUANTILES = (.05, .25, .50, .75, .95, .99)
METRICS = ("median_prediction", "median_signed_error", "median_absolute_error",
           "relative_median_absolute_error", "relative_underprediction",
           "capped_mae", "capped_rmse", "median_sample_eigenvalue_error",
           "plugin_tolerance_margin", "risk")


def mean_interval(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("need at least two finite independent observations")
    risk = float(x.mean())
    se = float(x.std(ddof=1)/np.sqrt(len(x)))
    half = t.ppf(.975, len(x)-1)*se
    return risk, se, risk-half, risk+half


def relative_frequency_interval(loss, raw):
    """Risk-adaptive event P(raw < .5 R): a conservative joint-band CI.

    Allocate 2.5% failure to a 97.5% Student risk CI and 2.5% to the uniform
    Dvoretzky-Kiefer-Wolfowitz CDF band. The union bound does not require
    independence of losses/predictions. Coverage inherits the Student interval's
    approximation; the CDF component is finite-sample distribution-free.
    Unlike a percentile bootstrap, the bound stays nondegenerate at zero events.
    """
    loss, raw = np.asarray(loss), np.asarray(raw)
    risk, se, _, _ = mean_interval(loss)
    if raw.shape != loss.shape or not np.isfinite(raw).all():
        raise ValueError('need finite matched loss and prediction vectors')
    half = t.ppf(1-.025/2, len(loss)-1)*se
    lower_risk, upper_risk = max(0.,risk-half), min(1.,risk+half)
    epsilon = np.sqrt(np.log(2/.025)/(2*len(raw)))
    lo = max(0.,float(np.mean(raw < .5*lower_risk))-epsilon)
    hi = min(1.,float(np.mean(raw < .5*upper_risk))+epsilon)
    return lo, hi


def wilson_interval(successes, count):
    if not 0 <= successes <= count or count < 1:
        raise ValueError("invalid binomial counts")
    p, z = successes/count, norm.ppf(.975)
    denominator = 1+z*z/count
    center = (p+z*z/(2*count))/denominator
    half = z*np.sqrt(p*(1-p)/count+z*z/(4*count**2))/denominator
    return max(0., center-half), min(1., center+half)


def quantile_interval(values, probability):
    """Conservative binomial order-statistic CI, valid also with capped ties.

    Open endpoints are +/-infinity when sample size cannot support a finite
    distribution-free bound. No percentile or rank jitter is introduced.
    """
    x = np.sort(np.asarray(values, dtype=float))
    if len(x) < 1 or not np.isfinite(x).all() or not 0 < probability < 1:
        raise ValueError("invalid quantile data or probability")
    # K=#observations below the population q-quantile; choose l,u such that
    # P(K<l)<=.025 and P(K>=u)<=.025. X_(0)=-inf; X_(M+1)=inf.
    lower_rank = int(binom.ppf(.025, len(x), probability))
    upper_rank = int(binom.ppf(.975, len(x), probability))+1
    lo = x[lower_rank-1] if lower_rank >= 1 else -np.inf
    hi = x[upper_rank-1] if upper_rank <= len(x) else np.inf
    return float(lo), float(hi)


def tolerance(risk):
    return np.maximum(.002, .20*np.asarray(risk))


def oracle_label(benchmark, lower, upper):
    """Check the entire risk CI, including both piecewise-linear breakpoints."""
    points = [lower, upper]
    points += [p for p in (.002/.20, benchmark) if lower <= p <= upper]
    points = np.asarray(points)
    margins = np.abs(benchmark-points)-tolerance(points)
    if np.max(margins) <= 0:
        return "within"
    if np.min(margins) > 0:
        return "outside"
    return "uncertain"


def joint_metrics(records, benchmark):
    """records[..., replication, (loss,raw,capped,gap)] -> paired statistics."""
    records = np.asarray(records, dtype=float)
    loss, raw, cap = records[..., 0], records[..., 1], records[..., 2]
    risk = np.mean(loss, axis=-1)
    if np.any(risk <= 0) or not np.isfinite(records).all():
        raise ValueError("positive risk and complete finite records required")
    absolute = np.abs(raw-risk[..., None])
    median = np.median(raw, axis=-1)
    mad = np.median(absolute, axis=-1)
    capped_error = cap-risk[..., None]
    return np.stack((median, median-risk, mad, mad/risk,
                     np.mean(raw < .5*risk[..., None], axis=-1),
                     np.mean(np.abs(capped_error), axis=-1),
                     np.sqrt(np.mean(capped_error**2, axis=-1)),
                     np.median(raw-benchmark, axis=-1),
                     mad-tolerance(risk), risk), axis=-1)


def bootstrap_metrics(records, benchmark, rng, repetitions, batch=100):
    """Resample rows, never columns. Recompute the risk in every bootstrap."""
    records = np.asarray(records)
    if records.ndim != 2 or records.shape[1] != 4 or len(records) < 2:
        raise ValueError("records must have shape (M>=2,4)")
    if repetitions < 100 or batch < 1:
        raise ValueError("need at least 100 bootstrap resamples and positive batch")
    output = []
    for start in range(0, repetitions, batch):
        indices = rng.integers(0, len(records), (min(batch, repetitions-start), len(records)))
        output.append(joint_metrics(records[indices], benchmark))
    return np.concatenate(output)


def percentile_interval(values):
    return tuple(np.quantile(values, [.025, .975], axis=0))


def calibration_edges(predictions):
    """Merge repeated quantile edges; never split observations with tied caps."""
    edges = np.unique(np.quantile(predictions, np.linspace(0, 1, 6)))
    if len(edges) == 1:
        return np.array([-np.inf, np.inf])
    return np.r_[-np.inf, edges[1:-1], np.inf]


def calibration_groups(predictions, edges):
    return np.searchsorted(edges[1:-1], predictions, side="left")
