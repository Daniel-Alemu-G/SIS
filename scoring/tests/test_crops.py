import unittest
import warnings

import numpy as np

from scoring.crops import best_crop, choose_season, hybrid_gate, selected_values, viability


class CropTests(unittest.TestCase):
    def test_hybrid_gate_hand_computation(self):
        scores = [0.9, 0.8, 0.95, 0.7, 0.9, 1]
        expected = np.sqrt(min(scores) * np.prod(np.asarray(scores, dtype=np.float32).astype(float)) ** (1 / 6))
        self.assertAlmostEqual(float(hybrid_gate(scores)), expected, places=7)

    def test_every_gate_can_force_zero(self):
        for index in range(6):
            scores = [1.0] * 6
            scores[index] = 0
            self.assertEqual(viability(scores, 1, slope_score=1, floor=0.8)[0], 0)

    def test_slope_forces_zero(self):
        self.assertEqual(viability([1, 1], 0.8, slope_score=0, floor=0.8)[0], 0)

    def test_unknown_gate_takes_precedence_over_zero(self):
        self.assertTrue(np.isnan(viability([np.nan, 0], 0.5, slope_score=0, floor=0.8)[0]))

    def test_missing_additional_is_not_zero(self):
        self.assertTrue(np.isnan(viability([1], np.nan, slope_score=1, floor=0.8)[0]))

    def test_perfect_additional_equals_gate(self):
        result, gate, _ = viability([0.4, 0.8], 1, slope_score=1, floor=0.8)
        self.assertEqual(result, gate)

    def test_fertility_floor_and_monotonicity(self):
        last = -1
        for score in np.linspace(0, 1, 11):
            result, gate, _ = viability([0.4, 0.8], score, slope_score=1, floor=0.8)
            self.assertGreaterEqual(float(result) + 1e-7, 0.8 * gate)
            self.assertGreaterEqual(result, last)
            last = result

    def test_log_geomean_does_not_underflow_for_small_positive_scores(self):
        gate = hybrid_gate([1e-7] * 10)
        self.assertGreater(gate, 0)
        self.assertAlmostEqual(float(gate), 1e-7, places=12)

    def test_season_best_score_then_higher_mean_then_row_order(self):
        choice = choose_season([[0.7, 1, 1], [0.8, 1, 1]], [[400, 500, 500], [300, 700, 500]])
        np.testing.assert_array_equal(choice, [1, 1, 0])

    def test_missing_comparison_does_not_silently_drop_a_season(self):
        self.assertEqual(choose_season([0.7, np.nan], [400, np.nan]), -1)

    def test_selected_temperature_uses_same_season(self):
        np.testing.assert_array_equal(selected_values([[10, 10], [25, 25]], np.array([0, 1])), [10, 25])

    def test_unknown_season_returns_nan(self):
        self.assertTrue(np.isnan(selected_values([10, 25], np.array(-1))))


class BestCropTests(unittest.TestCase):
    def test_best_crop_and_margin(self):
        result = best_crop({"teff": np.array([0.8]), "wheat": np.array([0.5])}, {"teff": 1, "wheat": 2}, minimum_score=0.25)
        self.assertEqual(result["best_crop"][0], 1)
        self.assertAlmostEqual(result["best_crop_margin"][0], 0.3, places=7)

    def test_below_minimum_is_no_crop_not_nodata(self):
        result = best_crop({"teff": 0.2}, {"teff": 1}, minimum_score=0.25)
        self.assertEqual(result["best_crop"], 0)
        self.assertAlmostEqual(float(result["best_crop_score"]), 0.2, places=7)

    def test_exact_minimum_remains_selected(self):
        self.assertEqual(best_crop({"teff": 0.25}, {"teff": 1}, minimum_score=0.25)["best_crop"], 1)

    def test_missing_crop_ignored_and_original_index_kept(self):
        result = best_crop({"teff": np.nan, "wheat": 0.8}, {"teff": 1, "wheat": 4}, minimum_score=0.25)
        self.assertEqual(result["best_crop"], 4)
        self.assertTrue(np.isnan(result["best_crop_margin"]))

    def test_all_missing_is_logical_255_and_null_scores(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            result = best_crop({"teff": np.nan, "wheat": np.nan}, {"teff": 1, "wheat": 2}, minimum_score=0.25)
        self.assertEqual(result["best_crop"], 255)
        self.assertTrue(np.isnan(result["best_crop_score"]))
        self.assertTrue(np.isnan(result["best_crop_margin"]))

    def test_exact_crop_tie_uses_original_index_not_dict_insertion(self):
        result = best_crop({"wheat": 0.8, "teff": 0.8}, {"teff": 1, "wheat": 2}, minimum_score=0.25)
        self.assertEqual(result["best_crop"], 1)
        self.assertEqual(result["best_crop_margin"], 0)

    def test_duplicate_indices_fail(self):
        with self.assertRaises(ValueError):
            best_crop({"a": 0.8, "b": 0.5}, {"a": 1, "b": 1}, minimum_score=0.25)
