import unittest
import warnings

import numpy as np

from scoring.soil import blend, calibrate_drought, data_weight, drought_score, limiting_factor, main_rainy_season, soil_quality

WEIGHTS = {"mean_weight": 0.5, "min_weight": 0.5}


class SoilTests(unittest.TestCase):
    def test_mean_min_combination(self):
        self.assertAlmostEqual(float(soil_quality([0.8, 0.4], **WEIGHTS)), 0.5, places=7)

    def test_missing_basic_propagates_even_with_present_es(self):
        self.assertTrue(np.isnan(soil_quality([np.nan, 0.4], [0.2], **WEIGHTS)))

    def test_missing_extra_is_ignored_and_zero_is_real(self):
        result = soil_quality([0.8, 0.4], [np.nan, 0.0], **WEIGHTS)
        self.assertAlmostEqual(float(result), 0.2, places=7)

    def test_all_missing_extra_returns_nan_without_warnings(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertTrue(np.isnan(soil_quality([0.8], [np.nan, np.nan], **WEIGHTS)))

    def test_all_missing_basic_has_no_reduction_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertTrue(np.isnan(soil_quality([np.nan], **WEIGHTS)))

    def test_blend_limits_and_missing_fallback(self):
        np.testing.assert_array_equal(blend([0.6] * 5, [0.2, 0.2, 0.2, np.nan, 0.2], [0, 1, 0.5, 1, np.nan]), np.array([0.6, 0.2, 0.4, 0.6, 0.6], dtype=np.float32))

    def test_missing_basic_cannot_be_replaced_by_detailed(self):
        self.assertTrue(np.isnan(blend(np.nan, 0.5, 1)))

    def test_data_weight_scales_by_available_active_factors(self):
        result = data_weight([0.8, 0.8, np.nan, 0], [[1, 1, 1, 1], [np.nan, 0, 1, 1]])
        np.testing.assert_array_equal(result, np.array([0.4, 0.8, 0, 0], dtype=np.float32))

    def test_no_active_es_produces_zero_weight(self):
        np.testing.assert_array_equal(data_weight([1, 0.5], []), [0, 0])

    def test_bad_confidence_and_shapes_fail(self):
        with self.assertRaises(ValueError):
            data_weight(1.2, [0.5])
        with self.assertRaises(ValueError):
            blend([0.5], [0.5, 0.5], [1])

    def test_limiting_factor_tie_uses_stable_order(self):
        result = limiting_factor({"n": [0.3, np.nan], "p": [0.3, np.nan]}, {"n": 1, "p": 2})
        self.assertEqual(result[0], 1)
        self.assertTrue(np.isnan(result[1]))

    def test_limiting_factor_ignores_missing_extra(self):
        result = limiting_factor({"n": 0.3, "zn": np.nan}, {"n": 1, "zn": 2}, valid=True)
        self.assertEqual(result, 1)

    def test_invalid_weights_do_not_change_formula(self):
        with self.assertRaises(ValueError):
            soil_quality([0.5], mean_weight=0.2, min_weight=0.2)

    def test_scalar_and_batch_agree_exactly(self):
        basic = [np.array([0.43, 0.27], dtype=np.float32), np.array([0.82, 0.39], dtype=np.float32)]
        extra = [np.array([0.19, np.nan], dtype=np.float32)]
        array = soil_quality(basic, extra, **WEIGHTS)
        for i in range(2):
            scalar = soil_quality([v[i] for v in basic], [v[i] for v in extra], **WEIGHTS)
            if np.isnan(array[i]):
                self.assertTrue(np.isnan(scalar))
            else:
                self.assertEqual(array[i], scalar)


class DroughtTests(unittest.TestCase):
    def test_main_season_and_exact_mean_tie(self):
        season, cv = main_rainy_season([600, 200, 500], [300, 700, 500], [0.2, 0.3, 0.4], [0.8, 0.7, 0.9])
        np.testing.assert_array_equal(season, [1, 2, 1])
        np.testing.assert_array_equal(cv, np.array([0.2, 0.7, 0.4], dtype=np.float32))

    def test_missing_mean_or_selected_cv_propagates(self):
        season, cv = main_rainy_season([np.nan, 500], [300, 200], [0.3, np.nan], [0.2, 0.2])
        self.assertEqual(season[0], 0)
        self.assertTrue(np.isnan(cv).all())

    def test_calibration_uses_only_valid_inside_selected_cv(self):
        result = calibrate_drought([100, 100, 100, 100], [50, 50, 50, 50], [0.1, 0.2, 0.3, 99], [9, 9, 9, 9], [1, 1, 1, 0], low_percentile=0, high_percentile=100, method="linear")
        self.assertEqual(result["valid_calibration_cells"], 3)
        self.assertAlmostEqual(result["cv_low"], 0.1, places=7)
        self.assertAlmostEqual(result["cv_high"], 0.3, places=7)

    def test_constant_cv_cannot_invent_cutoffs(self):
        with self.assertRaisesRegex(ValueError, "cut-offs"):
            calibrate_drought([100, 100], [50, 50], [0.2, 0.2], [0.5, 0.5], [1, 1], low_percentile=25, high_percentile=90, method="linear")

    def test_drought_penalty_bound(self):
        risk, overall = drought_score([0.5, 0.5, 0.5], [0, 0.5, 1], cv_low=0, cv_high=1, weight=0.2)
        np.testing.assert_array_equal(risk, [0, 0.5, 1])
        np.testing.assert_array_equal(overall, np.array([0.5, 0.45, 0.4], dtype=np.float32))
