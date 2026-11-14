"""Replication records: errors stay visible and stop inferential analysis."""
import numpy as np
import pandas as pd

from spectral_stability.e3_math import (LOSS_TOLERANCE, measurements,
                                       sample_covariances)
from .design import eigenvalues, seed_sequence

VALUE_COLUMNS = ['loss', 'plugin_raw', 'plugin_capped', 'sample_gap']


def simulate_configuration(config, condition, stage, replications):
    n, k, d = condition['n'], config['k'], config['d']
    sigma = np.diag(eigenvalues(condition['delta'], condition['t']))
    population = np.eye(d)[:, :k]
    sequence = seed_sequence(config, stage, condition['config_index'], 0)
    rng = np.random.default_rng(sequence)
    # Separate sampling failures from measurement failures; retain every row.
    try:
        covariances = sample_covariances(sigma, n, replications, rng, config['sampling'])
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
        raise RuntimeError(f"sampling failed for {condition}: {error}") from error
    try:
        values = np.column_stack(measurements(covariances, population, n, k))
        errors = ['']*replications
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        # Diagnose independently; never drop a difficult observation.
        values = np.full((replications, 4), np.nan)
        errors = []
        for index, covariance in enumerate(covariances):
            try:
                values[index] = measurements(covariance, population, n, k)
                errors.append('')
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
                errors.append(f'{type(error).__name__}: {error}')
    frame = pd.DataFrame(values, columns=VALUE_COLUMNS)
    frame['numerical_error'] = errors
    frame['nonfinite'] = ~np.isfinite(values).all(axis=1)
    frame['loss_out_of_bounds'] = ((values[:, 0] < -LOSS_TOLERANCE)
                                    | (values[:, 0] > 1+LOSS_TOLERANCE))
    frame['invalid_gap'] = values[:, 3] <= 0
    frame['negative_prediction'] = values[:, 1] < 0
    frame['failure'] = (frame.nonfinite | frame.loss_out_of_bounds | frame.invalid_gap
                        | frame.negative_prediction | (frame.numerical_error != ''))
    for key, value in condition.items():
        frame[key] = value
    frame['replication'] = np.arange(replications)
    frame['stage'] = stage
    frame['seed_spawn_key'] = '.'.join(map(str, sequence.spawn_key))
    return frame
