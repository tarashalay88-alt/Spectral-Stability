"""Mathematical utilities for PCA spectral-stability Experiment E1.

All data matrices use rows for observations and columns for variables.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh


FloatArray = NDArray[np.float64]


def _symmetric_matrix(matrix: ArrayLike, *, name: str) -> FloatArray:
    """Return ``matrix`` as a finite, real, symmetric square array."""
    result = np.asarray(matrix, dtype=float)
    if result.ndim != 2 or result.shape[0] != result.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    if not np.allclose(result, result.T, rtol=1e-12, atol=1e-12):
        raise ValueError(f"{name} must be symmetric")
    return result


def _validate_k(k: int, d: int) -> None:
    if not isinstance(k, (int, np.integer)) or isinstance(k, (bool, np.bool_)):
        raise TypeError("k must be an integer")
    if not 1 <= k < d:
        raise ValueError("k must satisfy 1 <= k < d")


def spiked_covariance(Q: ArrayLike, k: int, a: float) -> FloatArray:
    """Construct the rank-``k`` spiked covariance ``Q diag(lambda) Q.T``.

    Parameters
    ----------
    Q:
        Orthogonal eigenvector matrix with shape ``(d, d)``. Its first ``k``
        columns are the population principal directions.
    k:
        Principal-subspace dimension, satisfying ``1 <= k < d``.
    a:
        Positive spike strength. The leading eigenvalues are ``1 + a`` and
        all remaining eigenvalues are one.

    Returns
    -------
    ndarray, shape (d, d)
        Symmetric positive-definite population covariance matrix.
    """
    basis = np.asarray(Q, dtype=float)
    if basis.ndim != 2 or basis.shape[0] != basis.shape[1]:
        raise ValueError("Q must be a square matrix")
    if not np.all(np.isfinite(basis)):
        raise ValueError("Q must contain only finite values")

    d = basis.shape[0]
    _validate_k(k, d)
    if not np.isfinite(a) or a <= 0:
        raise ValueError("a must be a positive finite number")
    if not np.allclose(basis.T @ basis, np.eye(d), rtol=1e-10, atol=1e-10):
        raise ValueError("Q must be orthogonal")

    eigenvalues = np.ones(d)
    eigenvalues[:k] += a
    covariance = (basis * eigenvalues) @ basis.T
    return (covariance + covariance.T) / 2.0


def gaussian_observations(
    covariance: ArrayLike,
    n: int,
    rng: np.random.Generator,
) -> FloatArray:
    """Generate ``n`` independent zero-mean Gaussian observations.

    ``covariance`` has shape ``(d, d)`` and the returned array has shape
    ``(n, d)``. Randomness is supplied exclusively by the explicit ``rng``.
    """
    sigma = _symmetric_matrix(covariance, name="covariance")
    if not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be positive")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")

    return rng.multivariate_normal(
        mean=np.zeros(sigma.shape[0]),
        cov=sigma,
        size=n,
        check_valid="raise",
    )


def centered_sample_covariance(observations: ArrayLike) -> FloatArray:
    """Compute the centered sample covariance using denominator ``n - 1``.

    Parameters
    ----------
    observations:
        Data matrix with shape ``(n, d)``, where ``n >= 2``.

    Returns
    -------
    ndarray, shape (d, d)
        ``(X - X_bar).T @ (X - X_bar) / (n - 1)``.
    """
    data = np.asarray(observations, dtype=float)
    if data.ndim != 2:
        raise ValueError("observations must have shape (n, d)")
    if data.shape[0] < 2:
        raise ValueError("at least two observations are required")
    if data.shape[1] < 1:
        raise ValueError("at least one variable is required")
    if not np.all(np.isfinite(data)):
        raise ValueError("observations must contain only finite values")

    centered = data - data.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered / (data.shape[0] - 1)
    return (covariance + covariance.T) / 2.0


def top_k_eigenvectors(covariance: ArrayLike, k: int) -> FloatArray:
    """Return orthonormal eigenvectors for the ``k`` largest eigenvalues.

    The symmetric input has shape ``(d, d)`` and the output has shape
    ``(d, k)``. Columns are ordered by decreasing eigenvalue. Their individual
    signs are arbitrary, as usual for eigenvectors.
    """
    sigma = _symmetric_matrix(covariance, name="covariance")
    d = sigma.shape[0]
    _validate_k(k, d)

    eigenvalues, eigenvectors = eigh(sigma, subset_by_index=(d - k, d - 1))
    order = np.argsort(eigenvalues)[::-1]
    return eigenvectors[:, order]


def normalized_projector_loss(
    estimated_basis: ArrayLike,
    population_basis: ArrayLike,
) -> float:
    """Compute ``||P_hat - P||_F^2 / (2k)`` from two bases.

    Both inputs must have shape ``(d, k)`` and orthonormal columns. The result
    is invariant to signs and rotations within either represented subspace and
    lies in ``[0, 1]`` up to floating-point rounding.
    """
    estimated = np.asarray(estimated_basis, dtype=float)
    population = np.asarray(population_basis, dtype=float)
    if estimated.ndim != 2 or population.ndim != 2:
        raise ValueError("both bases must be two-dimensional")
    if estimated.shape != population.shape:
        raise ValueError("both bases must have the same shape (d, k)")
    d, k = estimated.shape
    if not 1 <= k <= d:
        raise ValueError("bases must have shape (d, k) with 1 <= k <= d")
    if not np.all(np.isfinite(estimated)) or not np.all(np.isfinite(population)):
        raise ValueError("bases must contain only finite values")

    identity = np.eye(k)
    if not np.allclose(estimated.T @ estimated, identity, rtol=1e-10, atol=1e-10):
        raise ValueError("estimated_basis must have orthonormal columns")
    if not np.allclose(population.T @ population, identity, rtol=1e-10, atol=1e-10):
        raise ValueError("population_basis must have orthonormal columns")

    estimated_projector = estimated @ estimated.T
    population_projector = population @ population.T
    difference = estimated_projector - population_projector
    loss = float(np.sum(difference * difference) / (2.0 * k))
    return float(np.clip(loss, 0.0, 1.0))


def first_order_benchmark(covariance: ArrayLike, n: int, k: int) -> float:
    """Compute the general first-order PCA projector-loss benchmark.

    Eigenvalues of the symmetric ``(d, d)`` covariance are sorted in
    decreasing order. The leading and trailing spectra must be separated:
    ``lambda_k > lambda_{k+1}``. The sample size must satisfy ``n >= 2``.
    """
    sigma = _symmetric_matrix(covariance, name="covariance")
    d = sigma.shape[0]
    _validate_k(k, d)
    if not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be at least 2")

    eigenvalues = eigh(sigma, eigvals_only=True)[::-1]
    if eigenvalues[-1] < -1e-12:
        raise ValueError("covariance must be positive semidefinite")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    leading = eigenvalues[:k, None]
    trailing = eigenvalues[None, k:]
    gaps = leading - trailing
    if np.any(gaps <= 0):
        raise ValueError("the top-k eigenspace must have a positive eigengap")

    terms = leading * trailing / gaps**2
    return float(np.sum(terms) / (k * (n - 1)))


def spiked_benchmark(d: int, n: int, k: int, a: float) -> float:
    """Compute ``(d-k)(1+a) / ((n-1)a^2)`` for the spiked model."""
    if not isinstance(d, (int, np.integer)) or isinstance(d, (bool, np.bool_)):
        raise TypeError("d must be an integer")
    if d < 2:
        raise ValueError("d must be at least 2")
    _validate_k(k, d)
    if not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be at least 2")
    if not np.isfinite(a) or a <= 0:
        raise ValueError("a must be a positive finite number")

    return float((d - k) * (1.0 + a) / ((n - 1) * a**2))