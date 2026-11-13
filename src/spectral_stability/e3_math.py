"""Gaussian covariance simulation and unmodified first-order risk estimates.

The E3 spectrum belongs in the experiment; these routines allow general d,k,n.
No eigengap floor, eigenvalue repair, or loss clipping is performed here.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import wishart

from PCA_stability import centered_sample_covariance, gaussian_observations

LOSS_TOLERANCE = 1e-12


def check_dimensions(d, k, n):
    for name, value in (("d", d), ("k", k), ("n", n)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
    if not 1 <= k < d or n < 2:
        raise ValueError("require 1 <= k < d and n >= 2")


def spectrum_benchmark(eigenvalues, n, k):
    """Formula for descending nonnegative eigenvalues, including singular gaps.

    A zero gap returns the actual IEEE inf/nan, to be recorded by the caller.
    A negative eigenvalue is an error, including a negative rounding excursion;
    it is never silently replaced by zero. Batch dimensions are supported.
    """
    lam = np.asarray(eigenvalues, dtype=float)
    if lam.ndim < 1:
        raise ValueError("eigenvalues must have a spectral axis")
    check_dimensions(lam.shape[-1], k, n)
    if not np.isfinite(lam).all() or (lam < 0).any():
        raise ValueError("eigenvalues must be finite and nonnegative")
    if (np.diff(lam, axis=-1) > 0).any():
        raise ValueError("eigenvalues must be descending")
    leading, trailing = lam[..., :k, None], lam[..., None, k:]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return np.sum(leading * trailing / (leading - trailing)**2,
                      axis=(-2, -1)) / (k * (n - 1))


def capped_benchmark(raw):
    """The explicitly separate bounded estimator. NaNs remain NaNs."""
    return np.minimum(1., np.maximum(0., raw))


def projector_loss(estimated, population):
    """Unclipped Frobenius loss for orthonormal (d,k) bases, including batches."""
    estimated, population = np.asarray(estimated), np.asarray(population)
    if estimated.shape[-2:] != population.shape[-2:]:
        raise ValueError("bases must have matching (d,k) shapes")
    d, k = estimated.shape[-2:]
    check_dimensions(d, k, 2)
    for basis in (estimated, population):
        if not np.isfinite(basis).all() or not np.allclose(
                np.swapaxes(basis, -1, -2) @ basis, np.eye(k), atol=1e-12, rtol=1e-12):
            raise ValueError("basis must be finite and orthonormal")
    difference = (estimated @ np.swapaxes(estimated, -1, -2)
                  - population @ np.swapaxes(population, -1, -2))
    return np.sum(difference**2, axis=(-2, -1)) / (2*k)


def sample_covariances(sigma, n, m, rng, method="wishart"):
    """Exact law of the centered covariance with df=n-1.

    Bartlett Wishart sampling is used when df>=d. For df<d, fall back to
    direct centered observations (a valid singular sample covariance).
    """
    sigma = np.asarray(sigma, dtype=float)
    if sigma.ndim != 2 or sigma.shape[0] != sigma.shape[1]:
        raise ValueError("sigma must be square")
    d = sigma.shape[0]
    check_dimensions(d, 1, n)
    if not np.isfinite(sigma).all() or not np.allclose(sigma, sigma.T, atol=1e-12):
        raise ValueError("sigma must be finite and symmetric")
    np.linalg.cholesky(sigma)  # population must be strictly positive definite
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or m < 1:
        raise ValueError("m must be a positive integer")
    if method not in ("gaussian", "wishart"):
        raise ValueError("unknown sampling method")
    if method == "gaussian" or n-1 < d:
        return np.stack([centered_sample_covariance(gaussian_observations(sigma, n, rng))
                         for _ in range(m)])
    samples = wishart.rvs(df=n-1, scale=sigma, size=m, random_state=rng) / (n-1)
    return np.asarray(samples).reshape(m, d, d)


def measurements(covariance, population_basis, n, k):
    """One eigendecomposition; retain unaltered spectra and projector loss."""
    check_dimensions(np.shape(covariance)[-1], k, n)
    lam, vectors = np.linalg.eigh(covariance)
    lam, vectors = lam[..., ::-1], vectors[..., ::-1]
    raw = spectrum_benchmark(lam, n, k)
    loss = projector_loss(vectors[..., :k], population_basis)
    return loss, raw, capped_benchmark(raw), lam[..., k-1]-lam[..., k]
