import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import subspace_angles

from PCA_stability import (centered_sample_covariance, first_order_benchmark,
                           normalized_projector_loss, top_k_eigenvectors)
from spectral_stability.e3_math import (LOSS_TOLERANCE, capped_benchmark, measurements,
    projector_loss, sample_covariances, spectrum_benchmark)
from experiments.e3.design import (configurations, eigenvalues, read_config, seed_sequence)
from experiments.e3.simulation import simulate_configuration


class Projectors(unittest.TestCase):
    def test_general_dimensions_projectors_angles_and_invariances(self):
        rng = np.random.default_rng(919)
        for d,k in [(2,1),(6,2),(9,7),(12,3)]:
            population,_ = np.linalg.qr(rng.normal(size=(d,k)))
            estimated,_ = np.linalg.qr(rng.normal(size=(d,k)))
            rotation,_ = np.linalg.qr(rng.normal(size=(k,k)))
            for basis in (population,estimated,top_k_eigenvectors(np.diag(np.arange(d,0,-1)),k)):
                p = basis@basis.T
                assert_allclose(p,p.T,atol=1e-12)
                assert_allclose(p@p,p,atol=1e-12)
                self.assertEqual(np.linalg.matrix_rank(p),k)
                self.assertAlmostEqual(np.trace(p),k)
            loss = projector_loss(estimated,population)
            self.assertGreaterEqual(loss,-LOSS_TOLERANCE)
            self.assertLessEqual(loss,1+LOSS_TOLERANCE)
            self.assertAlmostEqual(float(projector_loss(population,population)),0.)
            assert_allclose(projector_loss(estimated@rotation,-population),loss,atol=1e-12)
            assert_allclose(np.mean(np.sin(subspace_angles(estimated,population))**2),loss,atol=1e-12)
            assert_allclose(normalized_projector_loss(estimated,population),loss,atol=1e-12)

    def test_known_angles_and_orthogonal_loss(self):
        p = np.eye(4)[:,:2]
        self.assertAlmostEqual(float(projector_loss(np.eye(4)[:,2:],p)),1.)
        theta = .42
        q = np.array([[np.cos(theta),0],[0,1],[np.sin(theta),0],[0,0]])
        self.assertAlmostEqual(float(projector_loss(q,p)),np.sin(theta)**2/2)


class SpectraAndBenchmarks(unittest.TestCase):
    def test_entire_frozen_grid(self):
        config = read_config()
        self.assertEqual(len(configurations(config)),84)
        for c in configurations(config):
            lam = eigenvalues(c['delta'],c['t'])
            self.assertEqual(len(lam),12)
            self.assertTrue((lam>0).all())
            self.assertTrue((np.diff(lam)<=0).all())
            self.assertAlmostEqual(lam[2]-lam[3],c['delta'])
            self.assertAlmostEqual(lam.sum(),23.)
            self.assertEqual(lam[0],6.)
            basis = top_k_eigenvectors(np.diag(lam),3)
            assert_allclose(basis@basis.T,np.diag([1]*3+[0]*9),atol=1e-12)
            assert_allclose(spectrum_benchmark(lam,c['n'],3),
                            first_order_benchmark(np.diag(lam),c['n'],3),rtol=1e-12)

    def test_hand_computed_examples_and_n_minus_one(self):
        # (3*1)/(3-1)^2 /(n-1=10) = 3/40.
        self.assertAlmostEqual(float(spectrum_benchmark([3.,1.],11,1)),3/40)
        # top=(4,2), bottom=(1): sum=4/9+2=22/9; /[2*10]=11/90.
        self.assertAlmostEqual(float(spectrum_benchmark([4.,2.,1.],11,2)),11/90)
        x = np.array([[1.,2.],[3.,0.],[5.,4.]])
        assert_allclose(centered_sample_covariance(x),[[4.,2.],[2.,4.]])
        assert_allclose(centered_sample_covariance(x+100),centered_sample_covariance(x))

    def test_raw_extreme_no_floor_no_cap(self):
        a = float(spectrum_benchmark([1+1e-7,1.,.1],50,1))
        b = float(spectrum_benchmark([1+2e-7,1.,.1],50,1))
        self.assertGreater(a,1e12)
        assert_allclose(a/b,4.,rtol=1e-6)
        self.assertEqual(float(capped_benchmark(a)),1.)
        self.assertTrue(np.isinf(spectrum_benchmark([1.,1.],50,1)))
        self.assertTrue(np.isnan(capped_benchmark(np.nan)))

    def test_measurements_general_shapes(self):
        for d,k,n in [(3,1,9),(7,4,20)]:
            covariance = np.diag(np.arange(d,0,-1,dtype=float))
            loss,raw,cap,gap = measurements(covariance,np.eye(d)[:,:k],n,k)
            self.assertEqual(loss,0.)
            self.assertEqual(gap,1.)
            self.assertEqual(cap,min(1.,raw))

    def test_invalid_arguments(self):
        for lam,n,k in [([2.,1.],1,1),([2.,-1.],10,1),([1.,2.],10,1),([2.,1.],10,2)]:
            with self.assertRaises(ValueError):
                spectrum_benchmark(lam,n,k)
        with self.assertRaises(TypeError):
            spectrum_benchmark([2.,1.],True,1)


class Sampling(unittest.TestCase):
    def test_direct_gaussian_and_wishart_moments_and_risk(self):
        # Four contrasting protocol configurations, 3000 per method each.
        # 6 Monte Carlo SE tolerances control chance failures across many checks.
        m = 3000
        for index,(n,delta,t) in enumerate([(50,.05,0.),(50,.05,.8),(500,1.,0.),(2000,1.,.8)]):
            sigma = np.diag(eigenvalues(delta,t))
            samples,risks = [],[]
            for method,j in [('gaussian',0),('wishart',1)]:
                rng = np.random.default_rng(np.random.SeedSequence(48151,spawn_key=(index,j)))
                cov = sample_covariances(sigma,n,m,rng,method)
                samples.append(cov)
                for a,b in [(0,0),(2,2),(11,11),(0,2),(3,11)]:
                    expected_var = (sigma[a,a]*sigma[b,b]+sigma[a,b]**2)/(n-1)
                    self.assertLess(abs(cov[:,a,b].mean()-sigma[a,b]),6*np.sqrt(expected_var/m))
                    x = cov[:,a,b]
                    variance = x.var(ddof=1)
                    variance_se = np.sqrt((np.mean((x-x.mean())**4)-variance**2)/m)
                    self.assertLess(abs(variance-expected_var),6*variance_se)
                loss = measurements(cov,np.eye(12)[:,:3],n,3)[0]
                risks.append((loss.mean(),loss.std(ddof=1)/np.sqrt(m)))
            difference = abs(risks[0][0]-risks[1][0])
            combined_se = np.hypot(risks[0][1],risks[1][1])
            self.assertLess(difference,6*combined_se)
            print(f'Sampling validation n={n} delta={delta} t={t}: '
                  f'Gaussian risk={risks[0][0]:.6f}, Wishart={risks[1][0]:.6f}, '
                  f'difference/combined SE={difference/combined_se:.3f}',flush=True)

    def test_general_rotated_covariance_moments(self):
        sigma = np.array([[3.,.7,0],[.7,2.,.2],[0,.2,1.]])
        for j,method in enumerate(['gaussian','wishart']):
            sample = sample_covariances(sigma,35,5000,np.random.default_rng(201+j),method)
            se = np.sqrt((np.outer(sigma.diagonal(),sigma.diagonal())+sigma**2)/34/5000)
            self.assertTrue((np.abs(sample.mean(axis=0)-sigma)<6*se).all())

    def test_singular_sample_covariance_falls_back_to_gaussian(self):
        a = sample_covariances(np.eye(5),3,4,np.random.default_rng(5),'wishart')
        b = sample_covariances(np.eye(5),3,4,np.random.default_rng(5),'gaussian')
        assert_allclose(a,b,atol=0,rtol=0)
        self.assertTrue(all(np.linalg.matrix_rank(x)==2 for x in a))

    def test_seed_streams_unique_and_reproducible(self):
        config = read_config()
        keys = [tuple(seed_sequence(config,stage,c['config_index'],purpose).spawn_key)
                for stage in ['smoke','pilot','confirmation_reserved']
                for c in configurations(config) for purpose in (0,1,2)]
        self.assertEqual(len(keys),len(set(keys)))
        condition = configurations(config)[0]
        a = simulate_configuration(config,condition,'smoke',10)
        b = simulate_configuration(config,condition,'smoke',10)
        self.assertTrue(a.equals(b))
