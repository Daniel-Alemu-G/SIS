import copy
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from config import grid, stack, stack_v2
from pipeline.phase2_build import diagnostics, raster
from pipeline.phase2_build.tests.fixtures import calibrated_fixture, source_fixture
from scoring.engine import ScoringEngine, decode
from scoring.soil import blend


class DiagnosticFunctionsTests(unittest.TestCase):
    def test_band_statistics_use_inside_cells_without_clipping_outliers(self):
        values = np.array([[-1, 0], [2, grid.NODATA]], dtype=np.float32)
        inside = np.array([[True, True], [False, True]])
        record = {"name": "raw_input", "logical_dtype": "float32", "valid_min": 0, "valid_max": 1}
        stats = diagnostics.band_statistics(values, inside, record)
        self.assertEqual(stats["min"], -1)
        self.assertEqual(stats["max"], 0)
        self.assertEqual(stats["median"], -0.5)
        self.assertEqual(stats["missing_inside_cells"], 1)
        self.assertEqual(stats["below_sanity_range_cells"], 1)
        self.assertEqual(sum(stats["histogram"]["counts"]), 2)

    def test_mask_histogram_includes_meaningful_outside_zero(self):
        values = np.array([[1, 0], [1, 0]], dtype=np.float32)
        stats = diagnostics.band_statistics(values, values == 1, {"name": "ethiopia_mask", "logical_dtype": "uint8"})
        self.assertEqual(stats["count"], 4)
        self.assertEqual(stats["valid_inside_cells"], 2)
        self.assertEqual(stats["histogram"], {"kind": "categorical", "values": [0, 1], "counts": [2, 2]})

    def test_empty_constant_and_score_histograms_are_explicit(self):
        record = {"name": "soil_quality", "logical_dtype": "float32", "formula_version": stack_v2.FORMULA_VERSION}
        for values, kind in (([grid.NODATA], "empty"), ([0, 0], "constant"), ([0, 0.5, 1], "ten equal-width diagnostic bins")):
            array = np.asarray(values, dtype=np.float32)
            stats = diagnostics.band_statistics(array, np.ones(array.shape, dtype=bool), record)
            self.assertEqual(stats["histogram"]["kind"], kind)
            self.assertEqual(sum(stats["histogram"]["counts"]), stats["count"])
            json.dumps(stats, allow_nan=False)

    def test_nonfinite_physical_statistics_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            diagnostics.band_statistics(np.array([np.nan], dtype=np.float32), np.array([True]), {"name": "test"})

    def test_float32_class_cutoffs_enter_higher_class_and_missing_is_255(self):
        cutoffs = [25, 50, 70, 85]
        below = np.nextafter(np.float32(0.25), np.float32(0))
        values = np.array([below, 0.25, 0.5, 0.7, 0.85, 1, np.nan, grid.NODATA], dtype=np.float32)
        np.testing.assert_array_equal(diagnostics.score_classes(values, cutoffs), [0, 1, 2, 3, 4, 4, 255, 255])
        with self.assertRaises(ValueError):
            diagnostics.score_classes([0.5], [50, 25])

    def test_class_transitions_separate_missing_and_outside(self):
        a = np.array([0.3, 0.3, np.nan, 0.9], dtype=np.float32)
        b = np.array([0.6, np.nan, 0.1, 0.1], dtype=np.float32)
        result = diagnostics.class_changes(a, b, np.array([True, True, True, False]), [25, 50, 70, 85])
        self.assertEqual(result["compared_cells"], 1)
        self.assertEqual(result["changed_class_cells"], 1)
        self.assertEqual(result["became_missing_cells"], 1)
        self.assertEqual(result["became_valid_cells"], 1)
        self.assertEqual(result["transition_counts"][1][2], 1)

    def test_spearman_common_masks_ties_and_constants(self):
        a, b = np.array([0, 0, 0.4, 0.8], dtype=np.float32), np.array([0.2, 0.1, 0.4, 0.9], dtype=np.float32)
        pairs = diagnostics.crop_correlations({"a": a, "b": b, "constant": np.zeros_like(a)}, np.ones(a.shape, dtype=bool))
        self.assertAlmostEqual(pairs[0]["spearman"], float(spearmanr(a, b).statistic))
        self.assertIsNone(pairs[1]["spearman"])
        json.dumps(pairs, allow_nan=False)

    def test_spearman_uses_pairwise_known_cells_and_ignores_physical_nodata(self):
        a = np.array([0.1, 0.2, 0.3, grid.NODATA], dtype=np.float32)
        b = np.array([0.1, grid.NODATA, 0.3, 0.4], dtype=np.float32)
        result = diagnostics.crop_correlations({"a": a, "b": b}, np.ones(4, dtype=bool))[0]
        self.assertEqual(result["paired_valid_cells"], 2)
        self.assertAlmostEqual(result["spearman"], 1)
        self.assertTrue(result["above_specification_0_95"])

    def test_limiting_shares_follow_lookup_and_do_not_count_missing(self):
        shares = diagnostics.limiting_shares(np.array([1, 1, 2, grid.NODATA, np.nan], dtype=np.float32), {"gate": 1, "additional": 2}, {"gate"})
        self.assertEqual(shares["valid_limiting_cells"], 3)
        self.assertAlmostEqual(shares["additional_limiting_share"], 1 / 3)
        with self.assertRaisesRegex(ValueError, "legal"):
            diagnostics.limiting_shares(np.array([255], dtype=np.float32), {"gate": 1})

    def test_edge_effect_is_separated_from_basic_spatial_variation(self):
        inside = np.ones((2, 2), dtype=bool)
        weight = np.array([[0, 0.01], [0, 0.02]], dtype=np.float32)
        basic = np.array([[0.2, 0.8], [0.1, 0.9]], dtype=np.float32)
        blended = blend(basic, np.ones_like(basic), weight)
        result = diagnostics.edge_continuity(inside, weight, weight, {"soil_quality": basic}, {"soil_quality": blended}, fertility_floor=0.8)
        stats = result["bands"]["soil_quality"]
        self.assertEqual(result["horizontal_and_vertical_coverage_boundary_pairs"], 2)
        self.assertEqual(stats["zero_weight_bitwise_mismatches"], 0)
        self.assertLess(stats["absolute_confidence_only_jump"]["max"], 0.003)
        self.assertGreater(stats["absolute_raw_spatial_jump"]["max"], 0.7)

    def test_edge_hard_step_or_zero_weight_disagreement_fails(self):
        inside = np.ones((1, 2), dtype=bool)
        weight = np.array([[0, 0.01]], dtype=np.float32)
        basic = np.full((1, 2), 0.5, dtype=np.float32)
        changed = np.array([[0.5, 0.9]], dtype=np.float32)
        with self.assertRaisesRegex(ValueError, "theoretical blend bound"):
            diagnostics.edge_continuity(inside, weight, weight, {"soil_quality": basic}, {"soil_quality": changed}, fertility_floor=0.8)
        changed = basic.copy()
        changed[0, 0] += 0.001
        with self.assertRaises(ValueError):
            diagnostics.edge_continuity(inside, weight, weight, {"soil_quality": basic}, {"soil_quality": changed}, fertility_floor=0.8)

    def test_empty_edge_reports_null_not_nan(self):
        values = np.full((2, 2), 0.5, dtype=np.float32)
        result = diagnostics.edge_continuity(np.ones((2, 2), dtype=bool), np.zeros_like(values), np.zeros_like(values),
                                             {"soil_quality": values}, {"soil_quality": values}, fertility_floor=0.8)
        self.assertIsNone(result["bands"]["soil_quality"]["absolute_confidence_only_jump"]["max"])
        json.dumps(result, allow_nan=False)


class ConfiguredDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="phase2-diagnostics-", dir="/tmp/opencode")
        self.root = Path(self.temporary.name)
        self.source = source_fixture()
        self.engine, self.manifest = calibrated_fixture(self.root, self.source, six_crops=False)

    def tearDown(self):
        self.source = None
        self.temporary.cleanup()

    def test_sensitivity_uses_only_owner_values_without_mutating_baseline(self):
        original = copy.deepcopy(self.engine.config.params)
        variants = diagnostics.sensitivity_engines(self.engine)
        self.assertEqual(set(variants), {"mean_min_0.3_0.7", "fertility_floor_0.7", "drought_weight_0"})
        self.assertEqual(variants["fertility_floor_0.7"].config.params["viability"]["fertility_floor"], 0.7)
        self.assertEqual(self.engine.config.params, original)

    def geography(self, engine=None):
        engine = engine or self.engine
        values = raster.read_named(self.source, stack.band_indices_by_name(self.source), stack.CORE_BAND_NAMES)
        evaluated = engine.evaluate(values)
        return diagnostics.geographic_checks(self.source, engine, values["ethiopia_mask"] == 1,
            {crop.crop: evaluated.bands[f"viability_{crop.crop}"] for crop in engine.config.crops}, decode(values["elevation_m"]),
            {season: decode(values[f"{season}_mean_mm"]) for season in ("kiremt", "belg")},
            {"geographic_checks": {"sites": {"Addis Ababa": {"common_grid_row": 1, "common_grid_col": 1},
                                            "Danakil/Dallol diagnostic site": {"common_grid_row": 2, "common_grid_col": 2}}}})

    def test_capital_exclusion_remains_null_and_nearby_selection_is_disclosed(self):
        report = self.geography()
        original = report["sites"]["Original Addis Ababa reference"]
        nearby = report["sites"]["Central highlands near Addis Ababa"]
        self.assertIsNone(original["scores"]["viability_teff"])
        self.assertEqual((nearby["row"], nearby["col"]), (2, 1))
        self.assertIn("no scores were optimized", nearby["selection"])
        self.assertIn("1 km", nearby["resolution_label"])
        json.dumps(report, allow_nan=False)

    def test_skipped_sorghum_does_not_require_an_absent_rain_rule(self):
        config = replace(self.engine.config, requirements={**self.engine.config.requirements, "sorghum": {}})
        engine = ScoringEngine(config, self.engine.calibration)
        report = self.geography(engine)
        self.assertNotIn("Below 1200 m with adequate rain", report["sites"])
        self.assertIsNone(report["checks"]["below_1200_sorghum_above_barley"])
