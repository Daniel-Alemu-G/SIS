import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from config.grid import NODATA
from scoring.engine import ScoringEngine
from scoring.functions import rainfall_reliability_normal
from scoring.tests.fixtures import config_fixture, values_fixture


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="scoring-engine-", dir="/tmp/opencode")
        self.config = config_fixture(Path(self.temporary.name))
        self.calibration = {"cv_low": 0.1, "cv_high": 0.6}
        self.engine = ScoringEngine(self.config, self.calibration)

    def tearDown(self):
        self.temporary.cleanup()

    def test_array_and_scalar_scores_agree_bit_for_bit(self):
        points = [values_fixture() for _ in range(5)]
        points[1]["es_confidence"] = 0
        points[1]["es_zn"] = NODATA
        points[2]["sg_soc"] = NODATA
        points[3]["landcover"] = 80
        points[4]["belg_mean_mm"] = 500
        arrays = {name: np.array([point[name] for point in points], dtype=np.float32) for name in points[0]}
        result = self.engine.evaluate(arrays)
        for index, point in enumerate(points):
            cell = self.engine.score_cell(point)
            for name, array in result.bands.items():
                expected = None if not np.isfinite(array[index]) or (name == "best_crop" and array[index] == 255) else float(array[index])
                self.assertEqual(cell["scores"][name], expected, (index, name))

    def test_normal_reliability_has_percent_values_and_correct_source_season(self):
        point = values_fixture()
        point["belg_mean_mm"] = 700
        report = self.engine.score_cell(point)
        crop = report["crops"]["teff"]
        self.assertEqual(crop["chosen_season"], "belg")
        factor = next(f for f in crop["factors"] if f["id"] == "rain_reliability")
        expected = rainfall_reliability_normal(np.float32(point["belg_cv"]), threshold_ratio=0.75)
        self.assertEqual(factor["value"], expected)
        self.assertEqual(factor["unit"], "%")
        temperature = next(f for f in crop["factors"] if f["id"] == "temperature_season_c")
        self.assertEqual(temperature["value"], point["tmean_belg_c"])

    def test_no_es_returns_basic_and_null_detailed(self):
        point = values_fixture() | {"es_confidence": 0, "es_zn": None}
        report = self.engine.score_cell(point)
        self.assertEqual(report["soil"]["tier"], "basic")
        self.assertEqual(report["scores"]["data_weight"], 0)
        self.assertIsNone(report["scores"]["soil_quality_detailed"])
        self.assertEqual(report["scores"]["soil_quality"], report["scores"]["soil_quality_basic"])

    def test_dominant_tier_limiting_factor(self):
        point = values_fixture() | {"es_zn": 0, "es_confidence": 1}
        report = self.engine.score_cell(point)
        self.assertEqual(report["soil"]["limiting_factor"], "zinc")
        report = self.engine.score_cell(point | {"es_confidence": 0.1})
        self.assertEqual(report["soil"]["limiting_factor"], "organic_carbon")

    def test_soil_missing_basic_propagates(self):
        report = self.engine.score_cell(values_fixture() | {"sg_soc": NODATA})
        for name in ("soil_quality_basic", "soil_quality_detailed", "soil_quality", "soil_score_overall", "soil_limiting_factor", "viability_teff"):
            self.assertIsNone(report["scores"][name], name)

    def test_exclusions_and_unknown_landcover_do_not_hide_soil_quality(self):
        for landcover in (50, 80, NODATA):
            report = self.engine.score_cell(values_fixture() | {"landcover": landcover})
            self.assertIsNone(report["scores"]["viability_teff"])
            self.assertIsNone(report["scores"]["viability_teff_limiting"])
            self.assertIsNone(report["scores"]["best_crop"])
            self.assertIsNotNone(report["scores"]["soil_quality"])

    def test_current_cropland_is_not_a_filter(self):
        for landcover in (10, 20, 30, 40, 60, 95, 100):
            self.assertIsNotNone(self.engine.score_cell(values_fixture() | {"landcover": landcover})["scores"]["viability_teff"])

    def test_outside_mask_all_scores_and_weight_are_null(self):
        report = self.engine.score_cell(values_fixture() | {"ethiopia_mask": 0})
        self.assertTrue(all(value is None for value in report["scores"].values()))
        self.assertIsNone(report["soil"]["tier"])

    def test_slope_hard_zero_is_limiting_not_a_gate(self):
        report = self.engine.score_cell(values_fixture() | {"slope_mean_deg": 30})
        self.assertEqual(report["scores"]["viability_teff"], 0)
        self.assertEqual(report["crops"]["teff"]["limiting_factor"], "slope_deg")
        factor = next(f for f in report["crops"]["teff"]["factors"] if f["id"] == "slope_deg")
        self.assertEqual(factor["group"], "additional")

    def test_additional_factor_can_limit_when_gate_is_perfect(self):
        report = self.engine.score_cell(values_fixture() | {"sg_soc": 0, "slope_mean_deg": 0})
        self.assertEqual(report["crops"]["teff"]["gate"], 1)
        self.assertEqual(report["crops"]["teff"]["limiting_factor"], "organic_carbon")

    def test_missing_named_input_fails_loudly(self):
        values = values_fixture()
        del values["sg_soc"]
        with self.assertRaisesRegex(ValueError, "missing named"):
            self.engine.evaluate(values)

    def test_json_never_contains_nan(self):
        report = self.engine.score_cell(values_fixture() | {"es_zn": np.nan, "sg_soc": NODATA})
        text = json.dumps(report, allow_nan=False)
        self.assertIn("null", text)

    def test_data_weight_accounts_for_individual_missing_nutrients(self):
        zinc = self.config.es_factors[0]
        extra = replace(zinc, factor="copper", detailed_band="es_cu")
        config = replace(self.config, soil_factors=self.config.soil_factors + (extra,))
        engine = ScoringEngine(config, self.calibration)
        report = engine.score_cell(values_fixture() | {"es_cu": np.nan})
        self.assertAlmostEqual(report["scores"]["data_weight"], 0.4, places=7)
        self.assertEqual(report["soil"]["tier"], "partial")

    def test_no_active_es_always_uses_basic(self):
        config = replace(self.config, soil_factors=self.config.basic_factors)
        report = ScoringEngine(config, self.calibration).score_cell(values_fixture())
        self.assertEqual(report["scores"]["data_weight"], 0)
        self.assertIsNone(report["scores"]["soil_quality_detailed"])

    def test_missing_crop_reliability_row_uses_cv_fallback(self):
        requirements = {name: {key: row for key, row in rows.items() if key != "rain_reliability"} for name, rows in self.config.requirements.items()}
        report = ScoringEngine(replace(self.config, requirements=requirements), self.calibration).score_cell(values_fixture())
        factor = next(f for f in report["crops"]["teff"]["factors"] if f["id"] == "rain_reliability")
        self.assertEqual(factor["unit"], "fraction")
        self.assertEqual(factor["shape"], "ramp_down")
