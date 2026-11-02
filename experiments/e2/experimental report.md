# E2 — Spectral Geometry Beyond the Eigengap

## 1. Research Question

**Is the boundary eigengap alone sufficient to characterize finite-sample principal-subspace drift from the population, or does the rest of the covariance spectrum also matter?**

Given a population covariance matrix $\Sigma\in\mathbb{R}^{d\times d}$ with eigenvalues $\lambda_1\geq\cdots\geq\lambda_d\geq 0$, the boundary eigengap for the leading $k$-dimensional subspace is

$$
\delta_k=\lambda_k-\lambda_{k+1}.
$$

Spectral theory establishes the importance of eigenseparation for subspace stability. E2 considers covariance matrices with the same leading principal subspace, boundary eigengap, and total variance, but different trailing eigenvalues.

## 2. Hypothesis

**H2.** Two different covariance matrices with identical principal subspaces and boundary eigengaps can exhibit different principal-subspace losses because first-order estimation risk also depends on cross-boundary eigenvalue pairs.

For our chosen family of covariance matrices, the directional prediction is

$$
R(0.8;n)>R(0;n),
$$

where $R(t;n)$ is the expected estimation loss at spectral parameter $t$ and sample size $n$.

The distinction from Davis–Kahan is important: an eigengap-based perturbation bound controls error, but does not necessarily characterize its expected magnitude under random covariance-estimation noise.

## 3. Evaluation Metrics

Let $V_k\in\mathbb{R}^{d\times k}$ contain orthonormal leading eigenvectors of the population covariance matrix, and let $\widehat V_k$ contain the corresponding leading eigenvectors of the sample covariance matrix. Define the population and estimated projectors by

$$
P_k=V_kV_k^\top,\qquad \widehat P_k=\widehat V_k\widehat V_k^\top.
$$

### 3.1. Normalized squared Frobenius projector loss

Our primary metric is

$$
L_k=\frac{\|\widehat P_k-P_k\|_F^2}{2k}
=\frac{1}{k}\sum_{i=1}^k\sin^2\theta_i,
$$

where $\theta_i$ is the $i$-th principal angle between the estimated and population subspaces.

The loss measures **subspace discrepancy**, not individual eigenvector instability. It is invariant to eigenvector sign changes and orthogonal changes of basis within the same subspace.

### 3.2. Expected estimation risk

Define

$$
R(t;n)=\mathbb{E}[L_k(t;n)].
$$

With $M$ Monte Carlo replications, its empirical estimator is

$$
\widehat R(t;n)=\frac{1}{M}\sum_{m=1}^{M}L_k^{(m)}(t;n).
$$

### 3.3. Risk contrast

The primary comparison is the risk difference between the two endpoint spectra:

$$
\Delta_R(n)=R(0.8;n)-R(0;n),
$$

estimated by

$$
\widehat\Delta_R(n)=\widehat R(0.8;n)-\widehat R(0;n).
$$

## 4. Experimental Setup

### 4.1. Population covariance construction

We consider the family

$$
\Sigma(t)=\operatorname{diag}\!\left(
6,4,3,2,
\underbrace{1+t,\ldots,1+t}_{4},
\underbrace{1-t,\ldots,1-t}_{4}
\right).
$$

The experimental dimensions are $d=12$ and $k=3$, with

$$
t\in\{0,0.2,0.4,0.6,0.8\}.
$$

The primary comparison is $t=0$ versus $t=0.8$; intermediate values provide additional evidence about dependence on spectral redistribution.

| Quantity | Treatment |
|---|---|
| Ambient dimension $d$ | Fixed at $12$ |
| Target dimension $k$ | Fixed at $3$ |
| Leading eigenvalues | Fixed: $6,4,3$ |
| Fourth eigenvalue | Fixed: $2$ |
| Boundary eigengap $\delta_3$ | Fixed: $1$ |
| Total variance $\operatorname{tr}(\Sigma)$ | Fixed: $23$ |
| Leading population subspace | Fixed: $\operatorname{span}(e_1,e_2,e_3)$ |
| Trailing eigenvalue distribution | Varied with $t$ |
| Sample size $n$ | Varied |

For every $t$ in the grid,

$$
\delta_3=\lambda_3-\lambda_4=3-2=1,
\qquad
\operatorname{tr}(\Sigma(t))=23.
$$

### 4.2. Sampling and covariance estimation

The specified sampling model is

$$
X_1,\ldots,X_n\overset{\mathrm{i.i.d.}}{\sim}\mathcal{N}(0,\Sigma(t)).
$$

The intended sample covariance estimator is

$$
\widehat\Sigma(t)=\frac{1}{n-1}\sum_{r=1}^{n}(X_r-\overline X)(X_r-\overline X)^\top.
$$

Sample eigendecomposition provides $\widehat V_3$, which is used to construct $\widehat P_3$ and calculate $L_3$.

### 4.3. Actual simulation configuration

The saved pilot contains **500 Monte Carlo replications** for each combination of

$$
n\in\{250,500,1000,2000\}
\quad\text{and}\quad
t\in\{0,0.2,0.4,0.6,0.8\}.
$$

The saved replication-level data contain one loss per sample size, spectral parameter, and replication index: **10,000 losses**, with no missing observations or duplicate keys. The output files record random seed `20261031`.

The contrasts are computed from paired replication-level losses. The precise random-number coupling and covariance-estimation implementation should be confirmed in the experiment code.

## 5. Theoretical Benchmarks

### 5.1. First-order risk prediction

Under the fixed-dimensional Gaussian covariance-estimation model, the leading-order benchmark for normalized projector loss is

$$
B(\Sigma,n,k)=\frac{1}{k(n-1)}
\sum_{\substack{i\leq k\\j>k}}
\frac{\lambda_i\lambda_j}{(\lambda_i-\lambda_j)^2}.
$$

For a fixed population covariance with $\delta_k>0$, the approximation is

$$
R(\Sigma;n)=B(\Sigma,n,k)+o(n^{-1})
\qquad(n\to\infty).
$$

The benchmark depends on all cross-boundary eigenvalue pairs, not only on $\delta_k$.

For E2, direct evaluation gives

$$
B(\Sigma(0),n,3)=\frac{6.741852}{n-1},
\qquad
B(\Sigma(0.8),n,3)=\frac{10.939936}{n-1}.
$$

Therefore,

$$
\Delta_B(n)=B(\Sigma(0.8),n,3)-B(\Sigma(0),n,3)
=\frac{4.198085}{n-1}.
$$

### 5.2. Predicted numerical risks

| $n$ | $B(\Sigma(0),n,3)$ | $B(\Sigma(0.8),n,3)$ | $\Delta_B(n)$ |
|---:|---:|---:|---:|
| 250 | 0.027076 | 0.043935 | 0.016860 |
| 500 | 0.013511 | 0.021924 | 0.008413 |
| 1000 | 0.006749 | 0.010951 | 0.004202 |
| 2000 | 0.003373 | 0.005473 | 0.002100 |

These values are analytical first-order predictions, **not observed simulation losses**.

## 6. Empirical Results

The results below are from the saved E2 pilot. The means, paired contrasts, and uncertainty estimates have been checked against the replication-level observations.

### 6.1. Primary result — endpoint risk comparison

At $n=500$, comparing $t=0.8$ with $t=0$:

| Quantity | Observed value |
|---|---:|
| Baseline risk $\widehat R(0;500)$ | 0.013880 |
| Modified risk $\widehat R(0.8;500)$ | 0.023281 |
| Paired risk difference $\widehat\Delta_R(500)$ | **+0.009401** |
| Relative increase in mean loss | **+67.7%** |

The positive observed contrast indicates greater subspace-estimation risk for the redistributed spectrum at $n=500$, despite the unchanged eigengap.

**All primary endpoint comparisons:**

| $n$ | Empirical $\widehat R(0;n)$ | Empirical $\widehat R(0.8;n)$ | Empirical $\widehat\Delta_R(n)$ |
|---:|---:|---:|---:|
| 250 | 0.028002 | 0.046160 | +0.018157 |
| 500 | 0.013880 | 0.023281 | +0.009401 |
| 1000 | 0.006571 | 0.010869 | +0.004299 |
| 2000 | 0.003456 | 0.005525 | +0.002068 |

### 6.2. Monte Carlo uncertainty

The following are 95% normal-approximation Monte Carlo confidence intervals computed from paired differences.

| $n$ | Paired MCSE | 95% interval for $\widehat\Delta_R$ |
|---:|---:|---:|
| 250 | 0.000705 | [0.016775, 0.019540] |
| 500 | 0.000334 | [0.008746, 0.010055] |
| 1000 | 0.000148 | [0.004009, 0.004589] |
| 2000 | 0.000078 | [0.001916, 0.002220] |

Every interval lies above zero, supporting a positive expected-risk contrast under the specified simulation model.

The paired standard error is calculated using

$$
D_m=L_3^{(m)}(0.8)-L_3^{(m)}(0),
\qquad
\widehat{\operatorname{SE}}(\widehat\Delta_R)=\frac{s_D}{\sqrt M},
$$

where $s_D$ is the sample standard deviation of the 500 paired differences.

### 6.3. Behavior across the full spectral parameter grid

The experiment also investigated the intermediate values $t=0.2,0.4,0.6$. The table shows observed mean normalized projector loss.

| $n$ | $t=0$ | $t=0.2$ | $t=0.4$ | $t=0.6$ | $t=0.8$ |
|---:|---:|---:|---:|---:|---:|
| 250 | 0.028002 | 0.028744 | 0.031744 | 0.037013 | 0.046160 |
| 500 | 0.013880 | 0.014397 | 0.015853 | 0.018570 | 0.023281 |
| 1000 | 0.006571 | 0.006788 | 0.007456 | 0.008713 | 0.010869 |
| 2000 | 0.003456 | 0.003551 | 0.003868 | 0.004474 | 0.005525 |

![E2 Monte Carlo mean normalized projector loss versus spectral redistribution parameter, for four sample sizes.](figures/e2_mean_loss_vs_t.png)

*Figure 1. Mean normalized projector loss across spectral configurations; error bars show Monte Carlo uncertainty.*

For every tested sample size, mean estimation loss increases monotonically across the five $t$-settings. All 16 saved contrasts between $t>0$ and $t=0$ have positive unadjusted 95% Monte Carlo intervals.

### 6.4. Empirical contrast versus theoretical benchmark

This comparison directly evaluates the first-order risk prediction.

| $n$ | Empirical $\widehat\Delta_R$ | Theoretical $\Delta_B$ | Empirical − theory |
|---:|---:|---:|---:|
| 250 | 0.018157 | 0.016860 | +0.001297 |
| 500 | 0.009401 | 0.008413 | +0.000988 |
| 1000 | 0.004299 | 0.004202 | +0.000096 |
| 2000 | 0.002068 | 0.002100 | −0.000032 |

![Observed paired risk contrasts relative to t=0 and theoretical first-order contrasts across spectral parameters and sample sizes.](figures/e2_contrast_vs_theory.png)

*Figure 2. Observed paired contrasts versus first-order theoretical contrasts.*

![Difference between observed and theoretical first-order contrasts as a function of spectral redistribution.](figures/e2_prediction_error.png)

*Figure 3. Prediction error $\widehat\Delta_R-\Delta_B$ for the first-order contrast.*

The theoretical benchmark correctly predicts the **direction** of the contrast at all tested sample sizes. The absolute empirical–theoretical discrepancy is smaller at the larger sample sizes.

However, we should not claim perfect quantitative agreement. At $n=500$, the empirical contrast exceeds the first-order benchmark by approximately **2.96 paired Monte Carlo standard errors**, indicating a noticeable discrepancy at that sample size. At $n=1000$ and $n=2000$, the discrepancies are much smaller relative to their Monte Carlo standard errors.

These results are consistent with the usefulness of the first-order approximation at larger sample sizes, without establishing an exact finite-sample relationship.

## 7. Discussion and Limitations

### 7.1. Interpretation of results

The observed contrasts are positive at all tested sample sizes. Although the population covariance matrices share the same leading subspace, boundary eigengap, and total variance, their finite-sample principal-subspace estimation risks differ.

The first-order benchmark predicts the direction of the contrast; its quantitative agreement with the simulation improves at larger sample sizes. Estimation loss also decreases as sample size grows.

These numerical findings do **not** invalidate Davis–Kahan or other eigengap-based bounds. They show that a fixed boundary eigengap does not uniquely determine expected estimation risk under the particular covariance-estimation model studied.

### 7.2. Scientific limitations

| Limitation | Implication |
|---|---|
| Gaussian sampling model | Does not characterize heavy-tailed return distributions. |
| Independent observations | Excludes temporal dependence and volatility clustering. |
| Known synthetic population covariance | Does not address uncertainty about the true covariance in financial data. |
| Fixed $d=12$, $k=3$ | Does not establish identical behavior in high-dimensional asymptotic regimes. |
| Selected spectral family | Does not prove monotonicity for arbitrary spectral changes. |
| First-order perturbation approximation | Does not give exact finite-sample risk. |
| 500 Monte Carlo replications | Leaves finite simulation uncertainty. |
| Pilot experiment | Requires separation from any later confirmatory or robustness analysis. |

One further distinction is essential: **the true population leading subspace does not move in E2**. The observed discrepancy comes from covariance estimation, not from an underlying regime-induced subspace change. This makes E2 useful preparation for the later equity application, where estimation noise and genuine changes in market covariance geometry may coexist.

### 7.3. Conclusion

The observed pilot supports the hypothesis that covariance spectra with the same leading subspaces and boundary eigengaps can exhibit different finite-sample principal-subspace estimation risks.

The differences are consistent in direction with the first-order benchmark, and their magnitudes decrease as sample size grows. These findings motivate studying spectral information beyond the boundary eigengap when analyzing the stability of empirically estimated principal subspaces.

---

**Supporting pilot data:** [replication-level losses](data/e2_pilot_replications.csv) · [mean-risk summary](data/e2_pilot_means.csv) · [paired-contrast summary](data/e2_pilot_contrasts.csv). Numerical results are from the pilot; code-level reproducibility details (sampling implementation and random-number coupling) should be verified against the experiment source before final publication.
