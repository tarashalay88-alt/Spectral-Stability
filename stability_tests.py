"""Focused unit tests for :mod:`pca_stability`."""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from PCA_stability import (
    centered_sample_covariance,
    first_order_benchmark,
    gaussian_observations,
    normalized_projector_loss,
    spiked_benchmark,
    spiked_covariance,
    top_k_eigenvectors,
)


class TestSpikedCovariance(unittest.TestCase):
    def setUp(self) -> None:
        raw = np.random.default_rng(7).normal(size=(6, 6))
        self.Q, _ = np.linalg.qr(raw)

    def test_eigenstructure_and_dimensions(self) -> None:
        sigma = spiked_covariance(self.Q, k=2, a=3.0)

        self.assertEqual(sigma.shape, (6, 6))
        assert_allclose(sigma, sigma.T, atol=1e-14)
        assert_allclose(sigma @ self.Q[:, :2], 4.0 * self.Q[:, :2])
        assert_allclose(sigma @ self.Q[:, 2:], self.Q[:, 2:])
        assert_allclose(np.linalg.eigvalsh(sigma), [1, 1, 1, 1, 4, 4])


class TestSamplingAndSampleCovariance(unittest.TestCase):
    def test_sampling_shape_and_explicit_generator_reproducibility(self) -> None:
        sigma = np.diag([2.0, 1.0, 0.5])
        x1 = gaussian_observations(sigma, 8, np.random.default_rng(123))
        x2 = gaussian_observations(sigma, 8, np.random.default_rng(123))

        self.assertEqual(x1.shape, (8, 3))
        assert_allclose(x1, x2)

    def test_centering_and_n_minus_one_identity(self) -> None:
        x = np.array([[1.0, 2.0], [3.0, 0.0], [5.0, 4.0]])
        centered = x - x.mean(axis=0)
        expected = centered.T @ centered / 2

        actual = centered_sample_covariance(x)

        self.assertEqual(actual.shape, (2, 2))
        assert_allclose(actual, expected)
        assert_allclose(actual, np.cov(x, rowvar=False, ddof=1))


class TestPrincipalSubspaces(unittest.TestCase):
    def test_top_eigenvectors_and_projector_properties(self) -> None:
        theta = 0.37
        rotation = np.array(
            [
                [np.cos(theta), -np.sin(theta), 0.0],
                [np.sin(theta), np.cos(theta), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        sigma = rotation @ np.diag([5.0, 3.0, 1.0]) @ rotation.T

        basis = top_k_eigenvectors(sigma, k=2)
        projector = basis @ basis.T
        expected_projector = rotation[:, :2] @ rotation[:, :2].T

        self.assertEqual(basis.shape, (3, 2))
        assert_allclose(basis.T @ basis, np.eye(2), atol=1e-14)
        assert_allclose(projector, projector.T, atol=1e-14)
        assert_allclose(projector @ projector, projector, atol=1e-14)
        assert_allclose(projector, expected_projector, atol=1e-14)

    def test_projector_loss_identities_and_basis_invariance(self) -> None:
        population = np.eye(4)[:, :2]
        within_subspace_rotation = np.array([[0.0, -1.0], [1.0, 0.0]])

        self.assertAlmostEqual(
            normalized_projector_loss(
                population @ within_subspace_rotation,
                population,
            ),
            0.0,
        )
        self.assertAlmostEqual(
            normalized_projector_loss(np.eye(4)[:, 2:], population),
            1.0,
        )

        theta = 0.41
        estimated = np.column_stack(
            ([np.cos(theta), 0.0, np.sin(theta), 0.0], [0.0, 1.0, 0.0, 0.0])
        )
        # One principal angle is theta and the other is zero.
        self.assertAlmostEqual(
            normalized_projector_loss(estimated, population),
            np.sin(theta) ** 2 / 2,
        )


class TestBenchmarks(unittest.TestCase):
    def test_general_formula_directly(self) -> None:
        eigenvalues = np.array([5.0, 3.0, 2.0, 1.0])
        n, k = 11, 2
        expected_sum = sum(
            eigenvalues[i] * eigenvalues[j]
            / (eigenvalues[i] - eigenvalues[j]) ** 2
            for i in range(k)
            for j in range(k, len(eigenvalues))
        )
        expected = expected_sum / (k * (n - 1))

        self.assertAlmostEqual(
            first_order_benchmark(np.diag(eigenvalues), n=n, k=k),
            expected,
        )

    def test_spiked_and_general_benchmarks_are_equivalent(self) -> None:
        d, n, k, a = 7, 101, 3, 2.5
        raw = np.random.default_rng(42).normal(size=(d, d))
        Q, _ = np.linalg.qr(raw)
        sigma = spiked_covariance(Q, k=k, a=a)

        assert_allclose(
            first_order_benchmark(sigma, n=n, k=k),
            spiked_benchmark(d=d, n=n, k=k, a=a),
            rtol=1e-13,
            atol=0.0,
        )

    def test_general_benchmark_requires_cross_eigengap(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive eigengap"):
            first_order_benchmark(np.eye(3), n=10, k=1)


if __name__ == "__main__":
    unittest.main()