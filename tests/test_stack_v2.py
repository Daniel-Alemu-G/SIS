import copy
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from osgeo import gdal

from config import stack, stack_v2
from pipeline.phase2_build.tests.fixtures import v1_manifest_fixture, v2_manifest_fixture

gdal.UseExceptions()


class ScoredManifestTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="stack-v2-contract-", dir="/tmp/opencode")
        self.manifest = v2_manifest_fixture(Path(self.temporary.name))

    def tearDown(self):
        self.temporary.cleanup()

    def test_full_six_crop_stack_has_60_bands(self):
        records = stack_v2.validate_manifest(self.manifest, v1_manifest=v1_manifest_fixture())
        self.assertEqual(len(records), 60)
        self.assertEqual(tuple(records)[:38], stack.CORE_BAND_NAMES)
        self.assertEqual(tuple(records)[38:45], stack_v2.SOIL_BAND_NAMES)
        self.assertEqual(tuple(records)[45:49], ("viability_teff", "viability_teff_limiting", "viability_wheat", "viability_wheat_limiting"))
        self.assertEqual(tuple(records)[-3:], stack_v2.BEST_CROP_BAND_NAMES)

    def test_v1_contract_is_not_weakened_to_accept_v2(self):
        with self.assertRaisesRegex(stack.StackContractError, "stack_version"):
            stack.validate_manifest(self.manifest)

    def test_wrong_version_or_boolean_version_fails(self):
        for version in (1, True, 2.0):
            manifest = copy.deepcopy(self.manifest)
            manifest["stack_version"] = version
            with self.assertRaisesRegex(stack.StackContractError, "stack_version"):
                stack_v2.validate_manifest(manifest)

    def test_prefix_records_cannot_change_even_if_prefix_hash_is_replaced(self):
        self.manifest["bands"][0]["description"] = "Changed original description"
        self.manifest["parent_stack"]["band_records_sha256"] = stack_v2.metadata_sha256(self.manifest["bands"][:38])
        with self.assertRaisesRegex(stack.StackContractError, "identical"):
            stack_v2.validate_manifest(self.manifest, v1_manifest=v1_manifest_fixture())

    def test_prefix_hash_detects_changes_without_external_manifest(self):
        self.manifest["bands"][0]["description"] = "Changed description"
        with self.assertRaisesRegex(stack.StackContractError, "records changed"):
            stack_v2.validate_manifest(self.manifest)

    def test_reordered_new_bands_fail(self):
        self.manifest["bands"][38], self.manifest["bands"][39] = self.manifest["bands"][39], self.manifest["bands"][38]
        with self.assertRaisesRegex(stack.StackContractError, "v2 scoring contract"):
            stack_v2.validate_manifest(self.manifest)

    def test_missing_score_band_fails(self):
        self.manifest["bands"].pop()
        with self.assertRaisesRegex(stack.StackContractError, "band count"):
            stack_v2.validate_manifest(self.manifest)

    def test_categorical_encoding_and_mode_are_required(self):
        for key, bad_value in (("overview_resample", "average"), ("physical_dtype", "uint8"), ("nodata", 255)):
            manifest = copy.deepcopy(self.manifest)
            next(row for row in manifest["bands"] if row["name"] == "best_crop")[key] = bad_value
            with self.subTest(key=key), self.assertRaisesRegex(stack.StackContractError, key):
                stack_v2.validate_manifest(manifest)

    def test_limiting_factor_lookup_cannot_be_relabelled(self):
        next(row for row in self.manifest["bands"] if row["name"] == "soil_limiting_factor")["lookup"]["1"] = "invented"
        with self.assertRaisesRegex(stack.StackContractError, "lookup"):
            stack_v2.validate_manifest(self.manifest)

    def test_best_crop_ids_follow_csv_order(self):
        record = next(row for row in self.manifest["bands"] if row["name"] == "best_crop")
        self.assertEqual(record["lookup"], {"0": "none", "1": "teff", "2": "wheat", "3": "maize", "4": "sorghum", "5": "barley", "6": "chickpea"})
        self.assertEqual(record["logical_nodata"], 255)
        self.assertNotIn(255, record["valid_values"])

    def test_tie_policy_cannot_change(self):
        self.manifest["scoring"]["ties"]["crop"] = "Random tie breaking"
        with self.assertRaisesRegex(stack.StackContractError, "tie/missing"):
            stack_v2.validate_manifest(self.manifest)

    def test_missing_calibration_or_wrong_percentile_fails(self):
        for key, value in (("cv_high", float("nan")), ("low_percentile", 1), ("valid_calibration_cells", 0)):
            manifest = copy.deepcopy(self.manifest)
            manifest["scoring"]["calibration"][key] = value
            with self.subTest(key=key), self.assertRaises(stack.StackContractError):
                stack_v2.validate_manifest(manifest)

    def test_frozen_provisional_row_cannot_enter_manifest(self):
        self.manifest["scoring"]["configuration"]["soil_factors"][-1]["source"] = "PROVISIONAL"
        with self.assertRaisesRegex(stack.StackContractError, "PROVISIONAL"):
            stack_v2.validate_manifest(self.manifest)

    def test_configuration_hashes_must_match_build_inputs(self):
        self.manifest["scoring"]["configuration"]["inputs"]["params"]["sha256"] = "a" * 64
        with self.assertRaisesRegex(stack.StackContractError, "input hash"):
            stack_v2.validate_manifest(self.manifest)

    def test_approved_parameters_and_reliability_interpretation_are_frozen(self):
        engine = stack_v2.frozen_engine(self.manifest)
        self.assertEqual(engine.config.params["rainfall_reliability"]["threshold_ratio"], 0.75)
        self.assertIn("NOT probability of meeting a crop", self.manifest["scoring"]["reliability_meaning"])


class ScoredRasterHeaderTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="stack-v2-header-", dir="/tmp/opencode")
        self.manifest = v2_manifest_fixture(Path(self.temporary.name))
        self.dataset = gdal.GetDriverByName("MEM").Create("", 1, 1, len(self.manifest["bands"]), gdal.GDT_Float32)
        for key, value in {"stack_version": "2", "parent_stack_sha256": self.manifest["parent_stack"]["raster_sha256"],
                           "scoring_configuration_sha256": stack_v2.metadata_sha256(self.manifest["scoring"]["configuration"]),
                           "formula_version": stack_v2.FORMULA_VERSION,
                           "cv_low": str(self.manifest["scoring"]["calibration"]["cv_low"]),
                           "cv_high": str(self.manifest["scoring"]["calibration"]["cv_high"])}.items():
            self.dataset.SetMetadataItem(key, value)
        for row in self.manifest["bands"]:
            band = self.dataset.GetRasterBand(row["index"])
            band.SetDescription(row["name"])
            band.SetNoDataValue(stack.grid.NODATA)
            band.SetUnitType(row["unit"])
            for key in ("unit", "source", "source_file", "logical_dtype", "overview_resample"):
                band.SetMetadataItem(key, row[key])
            if "lookup" in row:
                band.SetMetadataItem("lookup", stack_v2.canonical_json(row["lookup"]))

    def tearDown(self):
        self.dataset = None
        self.temporary.cleanup()

    def validate(self, *, require_cog=False):
        with mock.patch.object(stack_v2.grid, "assert_matches_grid") as grid_check:
            result = stack_v2.assert_stack_matches_manifest(self.dataset, self.manifest, require_cog=require_cog)
            grid_check.assert_called_once_with(self.dataset)
            return result

    def test_header_matches_frozen_manifest(self):
        self.assertEqual(len(self.validate()), 60)

    def test_raster_version_mismatch_fails(self):
        self.dataset.SetMetadataItem("stack_version", "1")
        with self.assertRaisesRegex(stack.StackContractError, "stack_version"):
            self.validate()

    def test_configuration_tag_detects_changed_parameters(self):
        self.dataset.SetMetadataItem("scoring_configuration_sha256", "a" * 64)
        with self.assertRaisesRegex(stack.StackContractError, "scoring_configuration_sha256"):
            self.validate()

    def test_lookup_tag_must_match(self):
        index = stack.band_indices_by_name(self.dataset)["viability_teff_limiting"]
        self.dataset.GetRasterBand(index).SetMetadataItem("lookup", "{}")
        with self.assertRaisesRegex(stack.StackContractError, "lookup"):
            self.validate()

    def test_unit_and_nodata_tags_are_checked(self):
        self.dataset.GetRasterBand(39).SetUnitType("percent")
        with self.assertRaisesRegex(stack.StackContractError, "raster unit"):
            self.validate()

    def test_non_cog_is_not_accepted(self):
        with self.assertRaisesRegex(stack.StackContractError, "predictor-3 COG"):
            self.validate(require_cog=True)
