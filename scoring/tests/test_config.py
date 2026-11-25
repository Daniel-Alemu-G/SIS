import copy
import json
import tempfile
import unittest
from pathlib import Path

from scoring.config import ConfigurationError, from_metadata, load_config, validate_params
from scoring.tests.fixtures import manifest_fixture, params_fixture, requirements_fixture, soil_rows_fixture, write_fixture


class ConfigTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="scoring-config-", dir="/tmp/opencode")
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def load(self, **kwargs):
        return load_config(**write_fixture(self.root, **kwargs), manifest=manifest_fixture())

    def test_valid_configuration_and_original_csv_alias(self):
        config = self.load()
        self.assertEqual(len(config.soil_factors), 4)
        self.assertEqual([c.index for c in config.crops], [1, 2])
        temperature = config.requirements["teff"]["temperature_season_c"]
        self.assertEqual(temperature.csv_factor, "temperature")
        self.assertEqual(temperature.input_band, "tmean_<season>_c")
        self.assertEqual(len(config.inputs), 4)
        self.assertEqual(config.skipped, ())

    def test_alias_must_be_explicit(self):
        params = params_fixture()
        params["crop_factor_aliases"] = {}
        with self.assertRaisesRegex(ConfigurationError, "unknown"):
            self.load(params=params)

    def test_duplicate_soil_factor_fails(self):
        rows = soil_rows_fixture()
        rows.append(rows[0].copy())
        with self.assertRaisesRegex(ConfigurationError, "duplicate"):
            self.load(soil=rows)

    def test_duplicate_crop_alias_pair_fails(self):
        rows = requirements_fixture()
        extra = next(row.copy() for row in rows if row["crop"] == "teff" and row["factor"] == "temperature")
        extra["factor"] = "temperature_season_c"
        rows.append(extra)
        with self.assertRaisesRegex(ConfigurationError, "duplicate"):
            self.load(requirements=rows)

    def test_missing_required_crop_gate_skips_crop_without_renumbering(self):
        rows = [r for r in requirements_fixture() if (r["crop"], r["factor"]) != ("teff", "elevation_m")]
        config = self.load(requirements=rows)
        self.assertEqual([(c.crop, c.index) for c in config.crops], [("wheat", 2)])
        self.assertTrue(any(r["kind"] == "crop" for r in config.skipped))

    def test_tbd_factor_is_reported(self):
        rows = soil_rows_fixture()
        rows[-1]["a"] = "TBD"
        config = self.load(soil=rows)
        self.assertEqual(len(config.es_factors), 0)
        self.assertEqual(config.skipped[0]["reason"], "TBD value")

    def test_provisional_source_cannot_be_overridden_by_activation(self):
        rows = soil_rows_fixture()
        rows[-1]["source"] = "PROVISIONAL"
        config = self.load(soil=rows)
        self.assertFalse(config.es_factors)
        self.assertIn("PROVISIONAL", config.skipped[0]["reason"])

    def test_confirmed_units_resolve_stale_notes_without_changing_values(self):
        rows = soil_rows_fixture()
        rows[-1]["notes"] = "CSV unit unconfirmed; skip in v1"
        config = self.load(soil=rows)
        self.assertEqual(config.es_factors[0].a, 0)
        self.assertEqual(config.es_factors[0].b, 2)

    def test_unconfirmed_units_are_reported_and_skipped(self):
        rows = soil_rows_fixture()
        rows[-1]["notes"] = "CSV unit unconfirmed"
        params = params_fixture()
        params["confirmed_units"] = {}
        config = self.load(soil=rows, params=params)
        self.assertFalse(config.es_factors)

    def test_unit_mismatch_cannot_rescale_nutrients(self):
        rows = soil_rows_fixture()
        rows[1]["unit"] = "%"
        config = self.load(soil=rows)
        self.assertIn("unit incompatible", config.skipped[0]["reason"])

    def test_wrong_band_binding_fails(self):
        rows = soil_rows_fixture()
        rows[0]["band"] = "es_n"
        with self.assertRaisesRegex(ConfigurationError, "sg_phh2o"):
            self.load(soil=rows)

    def test_percent_reliability_cannot_use_raw_cv(self):
        params = params_fixture()
        params["rainfall_reliability"]["method"] = "cv"
        with self.assertRaisesRegex(ConfigurationError, "reliability shape/unit"):
            self.load(params=params)

    def test_crop_unit_mismatch_fails(self):
        rows = requirements_fixture()
        rows[0]["unit"] = "%"
        with self.assertRaisesRegex(ConfigurationError, "mismatch"):
            self.load(requirements=rows)

    def test_invalid_breakpoint_order_fails(self):
        rows = soil_rows_fixture()
        rows[0]["b"] = 10
        with self.assertRaisesRegex(ConfigurationError, "require a < b"):
            self.load(soil=rows)

    def test_ramp_extra_breakpoints_are_not_ignored(self):
        rows = soil_rows_fixture()
        rows[-1]["c"] = 5
        with self.assertRaisesRegex(ConfigurationError, "empty c/d"):
            self.load(soil=rows)

    def test_missing_activation_policy_fails(self):
        params = params_fixture()
        params["active_soil_factors"].remove("zinc")
        with self.assertRaisesRegex(ConfigurationError, "policy is missing"):
            self.load(params=params)

    def test_missing_manifest_band_fails(self):
        manifest = manifest_fixture()
        manifest["bands"] = [r for r in manifest["bands"] if r["name"] != "es_confidence"]
        with self.assertRaisesRegex(ConfigurationError, "lacks"):
            load_config(**write_fixture(self.root), manifest=manifest)

    def test_frozen_configuration_round_trip_without_source_reads(self):
        original = self.load()
        frozen = json.loads(json.dumps(original.metadata()))
        restored = from_metadata(frozen, manifest_fixture())
        self.assertEqual(json.loads(json.dumps(restored.metadata())), frozen)
        for path in self.root.iterdir():
            path.unlink()
        self.assertEqual(from_metadata(frozen, manifest_fixture()).soil_factors, original.soil_factors)

    def test_frozen_reported_missing_rule_can_be_restored(self):
        rows = soil_rows_fixture()
        rows[-1]["a"] = "TBD"
        config = self.load(soil=rows)
        restored = from_metadata(json.loads(json.dumps(config.metadata())), manifest_fixture())
        self.assertFalse(restored.es_factors)
        self.assertEqual(restored.skipped, config.skipped)

    def test_frozen_units_cannot_be_relabelled(self):
        frozen = self.load().metadata()
        frozen["band_units"]["kiremt_cv"] = "%"
        with self.assertRaisesRegex(ConfigurationError, "units"):
            from_metadata(frozen, manifest_fixture())

    def test_frozen_provisional_factor_is_rejected(self):
        frozen = self.load().metadata()
        frozen["soil_factors"][-1]["source"] = "PROVISIONAL"
        with self.assertRaisesRegex(ConfigurationError, "PROVISIONAL"):
            from_metadata(frozen, manifest_fixture())

    def test_frozen_crop_order_cannot_be_renumbered(self):
        frozen = self.load().metadata()
        frozen["crops"][-1]["index"] = 1
        with self.assertRaisesRegex(ConfigurationError, "increasing"):
            from_metadata(frozen, manifest_fixture())

    def test_frozen_alias_cannot_change(self):
        frozen = self.load().metadata()
        frozen["requirements"]["teff"]["temperature_season_c"]["csv_factor"] = "temprature"
        with self.assertRaisesRegex(ConfigurationError, "alias"):
            from_metadata(frozen, manifest_fixture())

    def test_frozen_missing_soil_rule_must_be_reported(self):
        frozen = self.load().metadata()
        frozen["soil_factors"].pop()
        with self.assertRaisesRegex(ConfigurationError, "without reporting"):
            from_metadata(frozen, manifest_fixture())

    def test_frozen_bad_fingerprint_fails(self):
        frozen = self.load().metadata()
        frozen["inputs"]["params"]["sha256"] = "not a hash"
        with self.assertRaisesRegex(ConfigurationError, "fingerprint"):
            from_metadata(frozen, manifest_fixture())


class ParameterTests(unittest.TestCase):
    def test_parameter_fixture_passes(self):
        validate_params(params_fixture())

    def test_missing_owner_approval_fails(self):
        params = params_fixture()
        params["approval"]["approved"] = False
        with self.assertRaisesRegex(ConfigurationError, "approval"):
            validate_params(params)

    def test_invalid_parameters_are_rejected(self):
        for section, key, value in (("soil_quality", "min_weight", 0.8), ("viability", "fertility_floor", 2),
                                    ("drought", "weight", float("nan")), ("drought", "cv_low_percentile", 100),
                                    ("blend", "detailed_label_min_weight", True)):
            with self.subTest(section=section, key=key):
                params = params_fixture()
                params[section][key] = value
                with self.assertRaises(ConfigurationError):
                    validate_params(params)

    def test_unknown_parameter_is_not_ignored(self):
        params = params_fixture()
        params["viability"]["typo"] = 1
        with self.assertRaisesRegex(ConfigurationError, "unknown"):
            validate_params(params)

    def test_duplicate_activation_fails(self):
        params = params_fixture()
        params["active_soil_factors"].append("ph")
        with self.assertRaisesRegex(ConfigurationError, "unique"):
            validate_params(params)

    def test_score_classes_must_be_ordered(self):
        params = params_fixture()
        params["score_classes"] = [50, 25]
        with self.assertRaisesRegex(ConfigurationError, "increasing"):
            validate_params(params)

    def test_invalid_sensitivity_pairs_fail(self):
        params = params_fixture()
        params["sensitivity"]["mean_min_weights"] = [[0.1, 0.2]]
        with self.assertRaisesRegex(ConfigurationError, "sum"):
            validate_params(params)

    def test_numeric_yaml_strings_are_not_silently_accepted(self):
        for section, key in (("soil_quality", "mean_weight"), ("rainfall_reliability", "threshold_ratio")):
            params = params_fixture()
            params[section][key] = str(params[section][key])
            with self.assertRaisesRegex(ConfigurationError, "numeric YAML"):
                validate_params(params)
