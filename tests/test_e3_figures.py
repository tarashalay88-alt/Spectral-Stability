"""Regression tests for E3 diagnostic figure scaling."""
import unittest

import pandas as pd

from experiments.e3.figures import serious_underprediction_color_limits


class SeriousUnderpredictionFigureTests(unittest.TestCase):
    def test_ineligible_large_rates_do_not_change_color_range(self):
        summary = pd.DataFrame({
            'serious_eligibility': ['established', 'established',
                                    'below_threshold', 'uncertain'],
            'serious_underprediction': [0.0, 0.445, 1.0, 0.91],
        })
        self.assertEqual(serious_underprediction_color_limits(summary), (0.0, 0.5))

    def test_zero_rate_still_has_nonzero_scale(self):
        summary = pd.DataFrame({
            'serious_eligibility': ['established', 'below_threshold'],
            'serious_underprediction': [0.0, 1.0],
        })
        self.assertEqual(serious_underprediction_color_limits(summary), (0.0, 0.1))

    def test_no_eligible_cells_uses_safe_fallback(self):
        summary = pd.DataFrame({
            'serious_eligibility': ['below_threshold', 'uncertain'],
            'serious_underprediction': [1.0, 0.25],
        })
        self.assertEqual(serious_underprediction_color_limits(summary), (0.0, 1.0))


if __name__ == '__main__':
    unittest.main()
