"""Transparent single-dataset inspection of the E1 mathematical core. NOT Monte Carlo E1."""
import numpy as np
from scipy.linalg import eigh, subspace_angles
from PCA_stability import (spiked_covariance, gaussian_observations, centered_sample_covariance,
                           top_k_eigenvectors, normalized_projector_loss,
                           first_order_benchmark, spiked_benchmark)

np.set_printoptions(precision=4, suppress=True, linewidth=120)
d, k, a, n = 4, 2, 2., 40
q = np.eye(d)
c, s = np.cos(0.5), np.sin(0.5)
q[0,0],q[0,2],q[2,0],q[2,2] = c,-s,s,c
c, s = np.cos(0.35), np.sin(0.35)
q[1,1],q[1,3],q[3,1],q[3,3] = c,-s,s,c
sigma = spiked_covariance(q,k,a)
x = gaussian_observations(sigma,n,np.random.default_rng(20261029))
sigma_hat = centered_sample_covariance(x)
u = q[:,:k]
u_hat = top_k_eigenvectors(sigma_hat,k)
p = u @ u.T
p_hat = u_hat @ u_hat.T
loss = normalized_projector_loss(u_hat,u)
angles = np.degrees(subspace_angles(u,u_hat))
print('SMALL VISUALIZATION (single draw; NOT E1 experiment)')
print('parameters:', f'd={d}, k={k}, a={a}, n={n}, seed=20261029')
print('Q =\n',q)
print('Population Sigma =\n',sigma)
print('Estimated Sigma_hat =\n',sigma_hat)
print('True eigenvalues =',np.linalg.eigvalsh(sigma)[::-1])
print('Estimated eigenvalues =',np.linalg.eigvalsh(sigma_hat)[::-1])
print('True projector P =\n',p)
print('Estimated projector P_hat =\n',p_hat)
print('||Sigma_hat-Sigma||_2 =',np.linalg.norm(sigma_hat-sigma,2))
print('Boundary eigengap =',a)
print('Principal angles degrees =',angles)
print('Normalized projector loss =',loss)
print('Overlap identity loss =',1-np.linalg.norm(u.T@u_hat,'fro')**2/k)
print('Theoretical first-order B =',first_order_benchmark(sigma,n,k))
print('Closed-form B =',spiked_benchmark(d,n,k,a))
print('Warning: single-draw loss vs B is not a direct prediction comparison. B predicts mean loss.')

print('\nPROJECT-SCALE SINGLE DRAW (no Monte Carlo):')
d,k,a,n = 50,3,2.,500
q = np.eye(d)
sigma = spiked_covariance(q,k,a)
x = gaussian_observations(sigma,n,np.random.default_rng(20261030))
sigma_hat = centered_sample_covariance(x)
u = q[:,:k]
u_hat=top_k_eigenvectors(sigma_hat,k)
print(f'd={d}, k={k}, n={n}; Sigma is diag(3,3,3,1,...,1)')
print('Sigma[0:6,0:6] =\n',sigma[:6,:6])
print('Sigma_hat[0:6,0:6] =\n',sigma_hat[:6,:6])
print('Top 7 estimated eigenvalues =',np.linalg.eigvalsh(sigma_hat)[-7:][::-1])
print('Principal angles degrees =',np.degrees(subspace_angles(u,u_hat)))
print('Normalized projector loss =',normalized_projector_loss(u_hat,u))
print('First-order benchmark B =',spiked_benchmark(d,n,k,a))