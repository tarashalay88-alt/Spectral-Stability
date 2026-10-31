"""E2 matched-spectrum Monte Carlo study. Run from the repository root.

Primary preregistered comparison: t=0 versus t=0.8 at n=500.
This script uses the existing PCA_stability.py functions; it does not
modify them or interpret results as market-regime evidence.
"""

import argparse
import csv
from pathlib import Path

import numpy as np

from PCA_stability import (
    centered_sample_covariance,
    first_order_benchmark,
    normalized_projector_loss,
    top_k_eigenvectors,
)

D, K = 12, 3
N_VALUES = (250, 500, 1000, 2000)
T_VALUES = (0.0, 0.2, 0.4, 0.6, 0.8)
PRIMARY_N, PRIMARY_T = 500, 0.8
PHASES = {
    "pilot": {"M": 500, "seed": 20261031},
    "final": {"M": 2000, "seed": 20261102},
}


def eigenvalues(t: float) -> np.ndarray:
    """Ordered population eigenvalues for E2, shape (12,)."""
    if not 0 <= t < 1:
        raise ValueError("t must satisfy 0 <= t < 1")
    return np.array([6., 4., 3., 2.] + [1. + t] * 4 + [1. - t] * 4)


def mean_and_mcse(values: np.ndarray) -> tuple[float, float]:
    """Sample mean and Monte Carlo standard error of independent replicates."""
    return float(np.mean(values)), float(np.std(values, ddof=1) / np.sqrt(len(values)))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run(phase: str) -> None:
    M, seed = PHASES[phase]["M"], PHASES[phase]["seed"]
    out = Path("results/e2")
    out.mkdir(parents=True, exist_ok=True)

    spectra = [eigenvalues(t) for t in T_VALUES]
    V_k = np.eye(D)[:, :K]  # Same true top-k subspace for every t.

    # Critical controls: same gap, leading eigenvalues, trace, and effective rank.
    for lam in spectra:
        assert np.isclose(lam[K - 1] - lam[K], 1.0)
        assert np.isclose(lam[0], 6.0)
        assert np.isclose(lam.sum(), 23.0)
        assert np.isclose(lam.sum() / lam[0], 23.0 / 6.0)

    # Different independent streams for n; within each replication, the same
    # standard-normal Z is reused for all spectra (paired common random numbers).
    streams = np.random.SeedSequence(seed).spawn(len(N_VALUES))
    summaries, contrasts, raw = [], [], []

    for n, stream in zip(N_VALUES, streams):
        rng = np.random.default_rng(stream)
        losses = np.empty((M, len(T_VALUES)))
        B = np.array([
            first_order_benchmark(np.diag(lam), n, K)
            for lam in spectra
        ])

        for m in range(M):
            Z = rng.standard_normal((n, D))
            for j, (t, lam) in enumerate(zip(T_VALUES, spectra)):
                # X = Z Lambda^(1/2): exact N(0, diag(lam)) sampling.
                X = Z * np.sqrt(lam)
                Sigma_hat = centered_sample_covariance(X)
                V_hat_k = top_k_eigenvectors(Sigma_hat, K)
                losses[m, j] = normalized_projector_loss(V_hat_k, V_k)
                raw.append({"n": n, "replication": m, "t": t,
                            "loss": float(losses[m, j])})

        for j, t in enumerate(T_VALUES):
            R_hat, mcse = mean_and_mcse(losses[:, j])
            summaries.append({
                "n": n, "t": t, "M": M, "seed": seed,
                "mean_loss": R_hat, "mcse": mcse,
                "ci_lower": R_hat - 1.96 * mcse,
                "ci_upper": R_hat + 1.96 * mcse,
                "benchmark_B": float(B[j]),
                "R_minus_B": R_hat - float(B[j]),
            })

            if j > 0:
                # Paired contrast: use the SD of differences, NOT independent SEs.
                differences = losses[:, j] - losses[:, 0]
                delta_R, delta_mcse = mean_and_mcse(differences)
                delta_B = float(B[j] - B[0])
                contrasts.append({
                    "n": n, "t": t, "reference_t": 0.0,
                    "M": M, "seed": seed,
                    "delta_R": delta_R,
                    "delta_mcse": delta_mcse,
                    "ci_lower": delta_R - 1.96 * delta_mcse,
                    "ci_upper": delta_R + 1.96 * delta_mcse,
                    "delta_B": delta_B,
                    "delta_R_minus_delta_B": delta_R - delta_B,
                })

        target = next(r for r in contrasts if r["n"] == n and r["t"] == PRIMARY_T)
        print(f"n={n:4d} | paired ΔR={target['delta_R']:.6f} "
              f"± {1.96 * target['delta_mcse']:.6f} "
              f"| theory ΔB={target['delta_B']:.6f}")

    write_csv(out / f"{phase}_means.csv", summaries)
    write_csv(out / f"{phase}_contrasts.csv", contrasts)
    write_csv(out / f"{phase}_replications.csv", raw)
    print(f"\nPrimary comparison: t=0 vs t={PRIMARY_T}, n={PRIMARY_N}")
    print(f"Saved {phase} results under {out}/ (M={M}, seed={seed}).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=PHASES, default="pilot")
    args = parser.parse_args()
    run(args.phase)
