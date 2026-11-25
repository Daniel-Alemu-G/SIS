import unittest

import numpy as np
from scipy.stats import norm

from scoring.functions import rainfall_reliability_normal, ramp_down, ramp_up, trapezoid


class FactorFunctionTests(unittest.TestCase):
    def test_ramp_up_breakpoints_and_bounds(self):
        np.testing.assert_array_equal(ramp_up([-2, 0, 1, 2, 4], 0, 2), [0, 0, 0.5, 1, 1])

    def test_ramp_down_breakpoints_and_bounds(self):
        np.testing.assert_array_equal(ramp_down([-2, 0, 1, 2, 4], 0, 2), [1, 1, 0.5, 0, 0])

    def test_trapezoid_breakpoints_and_plateau(self):
        np.testing.assert_array_equal(trapezoid([-1, 0, 1, 2, 3, 4, 5, 6, 7], 0, 2, 4, 6), [0, 0, 0.5, 1, 1, 1, 0.5, 0, 0])

    def test_triangle_is_valid(self):
        np.testing.assert_array_equal(trapezoid([0, 1, 2], 0, 1, 1, 2), [0, 1, 0])

    def test_scalar_and_array_float32_agree(self):
        array = ramp_up([0.7], 0, 3)
        scalar = ramp_up(0.7, 0, 3)
        self.assertEqual(scalar.shape, ())
        self.assertEqual(scalar.dtype, np.float32)
        self.assertEqual(array[0], scalar)

    def test_missing_and_nonfinite_are_preserved_as_missing(self):
        for function, args in ((ramp_up, (0, 1)), (ramp_down, (0, 1)), (trapezoid, (0, 1, 2, 3))):
            self.assertTrue(np.isnan(function([np.nan, np.inf, -np.inf], *args)).all())

    def test_invalid_parameters_fail(self):
        for a, b in ((1, 1), (2, 1), (np.nan, 2), (0, np.inf)):
            with self.assertRaises(ValueError):
                ramp_up(1, a, b)
        with self.assertRaises(ValueError):
            trapezoid(1, 0, 2, 1, 3)


class ReliabilityTests(unittest.TestCase):
    def test_matches_owner_normal_formula(self):
        cvs = np.array([0.1, 0.2, 0.3, 0.4])
        expected = (norm.cdf((1 - 0.75) / cvs) * 100).astype(np.float32)
        np.testing.assert_array_equal(rainfall_reliability_normal(cvs, threshold_ratio=0.75), expected)

    def test_zero_cv_owner_convention(self):
        self.assertEqual(rainfall_reliability_normal(0, threshold_ratio=1), 100)
        self.assertEqual(rainfall_reliability_normal(0, threshold_ratio=1.2), 0)

    def test_invalid_cv_is_not_imputed(self):
        self.assertTrue(np.isnan(rainfall_reliability_normal([-1, np.nan, np.inf], threshold_ratio=0.75)).all())

    def test_reliability_declines_as_cv_increases(self):
        result = rainfall_reliability_normal([0, 0.1, 0.2, 0.5, 1], threshold_ratio=0.75)
        self.assertTrue(np.all(np.diff(result) < 0))

    def test_scalar_and_array_agree(self):
        self.assertEqual(rainfall_reliability_normal(0.3, threshold_ratio=0.75), rainfall_reliability_normal([0.3], threshold_ratio=0.75)[0])

    def test_bad_threshold_ratio_fails(self):
        for ratio in (-1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                rainfall_reliability_normal(0.2, threshold_ratio=ratio)
