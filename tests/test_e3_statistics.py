import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import binom

import spectral_stability.e3_statistics as statistics
from spectral_stability.e3_statistics import (METRICS, bootstrap_metrics,
    calibration_edges, calibration_groups, joint_metrics, mean_interval,
    oracle_label, percentile_interval, quantile_interval, relative_frequency_interval,
    wilson_interval)


class Intervals(unittest.TestCase):
    def test_mcse_and_student_interval_hand_values(self):
        risk,se,lo,hi = mean_interval([0.,.2,.4,.6])
        self.assertAlmostEqual(risk,.3)
        self.assertAlmostEqual(se,np.sqrt(1/60))
        self.assertAlmostEqual(lo,.3-3.182446305284263*np.sqrt(1/60),places=10)
        self.assertAlmostEqual(hi,.3+3.182446305284263*np.sqrt(1/60),places=10)

    def test_wilson_boundaries_and_known_values(self):
        assert_allclose(wilson_interval(0,100),[0.,.03699349820698568],atol=1e-14)
        assert_allclose(wilson_interval(100,100),[.9630065017930143,1.],atol=1e-14)
        lo,hi = wilson_interval(50,100)
        self.assertAlmostEqual(lo+hi,1.)
        self.assertTrue(lo<.5<hi)

    def test_order_statistic_coverage_exact_and_ties(self):
        m = 1000
        x = np.arange(1,m+1,dtype=float)
        for probability in (.05,.25,.5,.75,.95,.99):
            lo,hi = quantile_interval(x,probability)
            lower_rank,upper_rank = int(lo),int(hi)
            coverage = binom.cdf(upper_rank-1,m,probability)-binom.cdf(lower_rank-1,m,probability)
            self.assertGreaterEqual(coverage,.95-1e-12)
        self.assertEqual(quantile_interval(np.ones(100),.5),(1.,1.))
        self.assertTrue(np.isposinf(quantile_interval(np.arange(40),.99)[1]))

    def test_mean_interval_and_bootstrap_mean_empirical_coverage(self):
        rng = np.random.default_rng(83972)
        hits_t,hits_bootstrap = 0,0
        for _ in range(400):
            loss = rng.uniform(0,1,60)
            _,_,lo,hi = mean_interval(loss)
            hits_t += lo<=.5<=hi
            records = np.column_stack([loss,.2+.2*loss,.2+.2*loss,1+loss])
            draws = bootstrap_metrics(records,.4,rng,500)
            lo,hi = percentile_interval(draws[:,METRICS.index('risk')])
            hits_bootstrap += lo<=.5<=hi
        # Binomial sampling error at 400 trials plus bootstrap finite-sample error.
        self.assertTrue(.89 < hits_t/400 < .99)
        self.assertTrue(.89 < hits_bootstrap/400 < .99)

    def test_oracle_classification_checks_whole_interval(self):
        self.assertEqual(oracle_label(.5,.49,.51),'within')
        self.assertEqual(oracle_label(.5,.2,.3),'outside')
        self.assertEqual(oracle_label(.5,.2,.8),'uncertain')
        self.assertEqual(oracle_label(.009,.006,.012),'uncertain')
        self.assertEqual(oracle_label(.009,.008,.010),'within')

    def test_relative_frequency_joint_interval_zero_and_coverage(self):
        lo,hi = relative_frequency_interval(np.full(1000,.3),np.ones(1000))
        self.assertEqual(lo,0.)
        self.assertGreater(hi,0.)
        self.assertLess(hi,.05)
        rng = np.random.default_rng(2671)
        hits = 0
        for _ in range(400):
            loss = rng.uniform(0,1,80)
            raw = .1+.8*loss  # R=.5, true P(raw<.25)=.1875; matched dependence
            lo,hi = relative_frequency_interval(loss,raw)
            hits += lo <= .1875 <= hi
        self.assertGreaterEqual(hits/400,.95)


class BootstrapAndTargets(unittest.TestCase):
    def setUp(self):
        # Skewness ensures median(L) != mean(L), separating the two targets.
        self.loss = .05+.85*np.linspace(0,1,80)**2
        self.raw = 3*self.loss+.1
        self.records = np.column_stack([self.loss,self.raw,np.minimum(1,self.raw),5-self.loss])

    def test_bootstrap_preserves_complete_records_and_recomputes_risk(self):
        seen = []
        original = joint_metrics
        def tracked(records,benchmark):
            assert_allclose(records[...,1],3*records[...,0]+.1)
            assert_allclose(records[...,2],np.minimum(1,records[...,1]))
            assert_allclose(records[...,3],5-records[...,0])
            seen.append(records.copy())
            return original(records,benchmark)
        with patch.object(statistics,'joint_metrics',side_effect=tracked):
            draws = bootstrap_metrics(self.records,.5,np.random.default_rng(234),200,batch=40)
        sampled = np.concatenate(seen)
        for index in (0,50,199):
            rows = sampled[index]
            risk = rows[:,0].mean()
            self.assertAlmostEqual(draws[index,METRICS.index('risk')],risk)
            self.assertAlmostEqual(draws[index,METRICS.index('median_absolute_error')],
                                   np.median(np.abs(rows[:,1]-risk)))
            self.assertAlmostEqual(draws[index,METRICS.index('relative_underprediction')],
                                   np.mean(rows[:,1]<.5*risk))
        self.assertGreater(np.std(draws[:,-1]),0)

    def test_bootstrap_risk_sd_matches_mcse(self):
        draws = bootstrap_metrics(self.records,.5,np.random.default_rng(811),4000)
        _,se,_,_ = mean_interval(self.loss)
        self.assertLess(abs(draws[:,-1].std(ddof=1)/se-1),.06)

    def test_target_is_population_risk_not_realized_loss(self):
        metric = joint_metrics(self.records,.5)
        risk = self.loss.mean()
        self.assertAlmostEqual(metric[METRICS.index('median_absolute_error')],
                               np.median(abs(self.raw-risk)))
        self.assertNotAlmostEqual(metric[METRICS.index('median_absolute_error')],
                                  np.median(abs(self.raw-self.loss)))
        self.assertAlmostEqual(metric[METRICS.index('capped_mae')],
                               np.abs(self.records[:,2]-risk).mean())
        self.assertTrue((np.abs(self.records[:,2]-risk)<=np.abs(self.raw-risk)+1e-14).all())
        self.assertAlmostEqual(metric[METRICS.index('median_signed_error')],
            metric[METRICS.index('median_sample_eigenvalue_error')]+.5-risk)

    def test_invalid_records_refused_not_dropped(self):
        broken = self.records.copy()
        broken[0,1] = np.inf
        with self.assertRaises(ValueError):
            bootstrap_metrics(broken,.5,np.random.default_rng(3),200)

    def test_tied_caps_merged_and_never_split(self):
        predictions = np.r_[np.linspace(.01,.8,200),np.ones(800)]
        edges = calibration_edges(predictions)
        groups = calibration_groups(predictions,edges)
        self.assertLess(len(edges)-1,5)
        self.assertEqual(len(np.unique(groups[predictions==1])),1)
        self.assertEqual(sum((groups==g).sum() for g in np.unique(groups)),1000)
        all_tied = calibration_edges(np.ones(100))
        self.assertEqual(len(all_tied)-1,1)
        self.assertTrue((calibration_groups(np.ones(100),all_tied)==0).all())


class FailureRecording(unittest.TestCase):
    def test_nonfinite_and_measurement_failures_are_retained(self):
        from experiments.e3.design import read_config,configurations
        from experiments.e3.simulation import simulate_configuration
        config = read_config()
        c = configurations(config)[0]
        with patch('experiments.e3.simulation.measurements',return_value=(
                np.array([.2,.3]),np.array([np.inf,np.nan]),np.array([1.,np.nan]),np.array([0.,.1]))):
            result = simulate_configuration(config,c,'smoke',2)
        self.assertEqual(len(result),2)
        self.assertTrue(result.failure.all())
        self.assertTrue(np.isinf(result.plugin_raw.iloc[0]))
        with patch('experiments.e3.simulation.measurements',side_effect=np.linalg.LinAlgError('test failure')):
            result = simulate_configuration(config,c,'smoke',2)
        self.assertEqual(len(result),2)
        self.assertTrue(result.failure.all())
        self.assertTrue(result.numerical_error.str.contains('test failure').all())

    def test_out_of_range_loss_flagged_without_clipping(self):
        from experiments.e3.design import read_config,configurations
        from experiments.e3.simulation import simulate_configuration
        config = read_config()
        with patch('experiments.e3.simulation.measurements',return_value=(
                np.array([1.01]),np.array([5.]),np.array([1.]),np.array([.1]))):
            result = simulate_configuration(config,configurations(config)[0],'smoke',1)
        self.assertEqual(result.loss.iloc[0],1.01)
        self.assertTrue(result.loss_out_of_bounds.iloc[0])
