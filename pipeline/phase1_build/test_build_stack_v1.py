"""Small synthetic regressions for lossless assembly and publication guards."""

import copy
import json
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest import mock

import numpy as np
from osgeo import gdal

import build_stack_v1 as build

gdal.UseExceptions()


def phosphorus_fixture():
    return {"status": "passed", "stack_eligible": True, "publication_ready": True, "publication_blockers": [],
            "source": {"upstream_product": "AfricanSoilGrids"},
            "bands": [{"name": "sg_p", "unit": "mg/100kg", "transform_applied": "none"}],
            "owner_approval": {"product": "AfricanSoilGrids", "unit": "mg/100kg"}}


class ApprovalTests(unittest.TestCase):
    def setUp(self):
        self.report = phosphorus_fixture()

    def test_approved_african_product_passes(self):
        build.validate_phosphorus_approval(self.report)

    def test_pending_status_fails(self):
        self.report["status"] = "staged_pending_provenance"
        with self.assertRaisesRegex(ValueError, "passed"):
            build.validate_phosphorus_approval(self.report)

    def test_false_stack_eligibility_fails(self):
        self.report["stack_eligible"] = False
        with self.assertRaisesRegex(ValueError, "inclusion approval"):
            build.validate_phosphorus_approval(self.report)

    def test_publication_blocker_fails(self):
        self.report["publication_blockers"] = ["Unverified source"]
        with self.assertRaisesRegex(ValueError, "publication blocker"):
            build.validate_phosphorus_approval(self.report)

    def test_not_publication_ready_fails(self):
        self.report["publication_ready"] = False
        with self.assertRaisesRegex(ValueError, "publication blocker"):
            build.validate_phosphorus_approval(self.report)

    def test_global_soilgrids_is_not_the_approved_source(self):
        self.report["source"]["upstream_product"] = "SoilGrids 2.0"
        with self.assertRaisesRegex(ValueError, "AfricanSoilGrids"):
            build.validate_phosphorus_approval(self.report)

    def test_mg_per_kg_cannot_replace_mg_per_100kg(self):
        self.report["bands"][0]["unit"] = "mg/kg"
        with self.assertRaisesRegex(ValueError, "as-is in mg/100kg"):
            build.validate_phosphorus_approval(self.report)

    def test_back_transform_cannot_be_invented(self):
        self.report["bands"][0]["transform_applied"] = "exp(value)"
        with self.assertRaisesRegex(ValueError, "without a transform"):
            build.validate_phosphorus_approval(self.report)

    def test_unmatched_owner_authority_fails(self):
        self.report["owner_approval"]["unit"] = "mg/kg"
        with self.assertRaisesRegex(ValueError, "provenance authority"):
            build.validate_phosphorus_approval(self.report)


class CopyTests(unittest.TestCase):
    def test_valid_values_and_signed_zero_are_bitwise_preserved(self):
        values = np.array([[0.0, -0.0, 174.62637, 1596.8274, 9091.994]], dtype=np.float32)
        result = build.masked_copy(values, np.ones_like(values, dtype=bool), build.grid.NODATA)
        self.assertEqual(build.bitwise_mismatches(values, result), 0)

    def test_nonfinite_values_and_tagged_missing_become_canonical_nodata(self):
        values = np.array([[np.nan, np.inf, -np.inf, build.grid.NODATA, 5]], dtype=np.float32)
        result = build.masked_copy(values, np.ones_like(values, dtype=bool), build.grid.NODATA)
        np.testing.assert_array_equal(result, [[build.grid.NODATA] * 4 + [5]])

    def test_outside_values_are_masked_without_changing_inside_outliers(self):
        values = np.array([[9091.994, 500]], dtype=np.float32)
        result = build.masked_copy(values, np.array([[True, False]]), build.grid.NODATA)
        self.assertEqual(result[0, 0].view(np.uint32), values[0, 0].view(np.uint32))
        self.assertEqual(result[0, 1], build.grid.NODATA)

    def test_zero_confidence_remains_valid_zero(self):
        values = np.array([[0, 0.5]], dtype=np.float32)
        np.testing.assert_array_equal(build.masked_copy(values, np.ones_like(values, dtype=bool), np.nan), values)

    def test_landcover_conversion_preserves_codes_and_missingness(self):
        values = np.array([[10, 95, 0]], dtype=np.uint8)
        result = build.masked_copy(values, np.ones_like(values, dtype=bool), 0)
        self.assertEqual(result.dtype, np.float32)
        np.testing.assert_array_equal(result, [[10, 95, build.grid.NODATA]])

    def test_shape_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, "same shape"):
            build.masked_copy(np.ones((1, 2)), np.ones((1, 1), dtype=bool), None)

    def test_nonboolean_mask_fails(self):
        with self.assertRaisesRegex(ValueError, "boolean"):
            build.masked_copy(np.ones((1, 1)), np.ones((1, 1)), None)

    def test_bitwise_comparison_rejects_dtype_mismatch(self):
        with self.assertRaisesRegex(ValueError, "float32"):
            build.bitwise_mismatches(np.ones((1, 1)), np.ones((1, 1), dtype=np.float32))


class LogicalTests(unittest.TestCase):
    def test_valid_mask_and_ids_pass(self):
        build.validate_logical_values("ethiopia_mask", np.array([[0, 1]], dtype=np.float32), set(), set())
        build.validate_logical_values("woreda_id", np.array([[0, 12, build.grid.NODATA]], dtype=np.float32), {12}, set())

    def test_mask_cannot_have_nodata(self):
        with self.assertRaisesRegex(ValueError, "meaningful 0/1"):
            build.validate_logical_values("ethiopia_mask", np.array([[build.grid.NODATA]], dtype=np.float32), set(), set())

    def test_fractional_ids_fail(self):
        with self.assertRaisesRegex(ValueError, "fractional"):
            build.validate_logical_values("woreda_id", np.array([[12.5]], dtype=np.float32), {12}, set())

    def test_unknown_ids_fail(self):
        with self.assertRaisesRegex(ValueError, "illegal"):
            build.validate_logical_values("woreda_id", np.array([[13]], dtype=np.float32), {12}, set())

    def test_unapproved_landcover_codes_fail(self):
        with self.assertRaisesRegex(ValueError, "illegal"):
            build.validate_logical_values("landcover", np.array([[11]], dtype=np.float32), set(), {10, 20})


class OverviewTests(unittest.TestCase):
    def fixture(self):
        names = ("elevation_m", "landcover", "woreda_id", "ethiopia_mask")
        dataset = gdal.GetDriverByName("MEM").Create("", 4, 4, 4, gdal.GDT_Float32)
        arrays = [
            np.array([[1, 1, 3, 3], [1, 2, 3, 4], [5, 5, 7, 7], [5, 6, 7, 8]], dtype=np.float32),
            np.array([[10, 10, 30, 30], [10, 20, 30, 40], [60, 60, 80, 80], [60, 95, 80, 90]], dtype=np.float32),
            np.array([[1, 1, 3, 3], [1, 2, 3, 4], [5, 5, 7, 7], [5, 6, 7, 8]], dtype=np.float32),
            np.array([[1, 1, 1, 1], [1, 0, 1, 0], [1, 1, 1, 1], [1, 0, 1, 0]], dtype=np.float32),
        ]
        for index, (name, values) in enumerate(zip(names, arrays), start=1):
            band = dataset.GetRasterBand(index)
            band.SetDescription(name)
            band.SetNoDataValue(build.grid.NODATA)
            band.WriteArray(values)
        return dataset

    def test_continuous_average_and_categorical_mode_are_distinct(self):
        dataset = self.fixture()
        indices = build.stack.band_indices_by_name(dataset)
        build.apply_overview_policy(dataset, indices, [2])
        np.testing.assert_array_equal(dataset.GetRasterBand(indices["elevation_m"]).GetOverview(0).ReadAsArray(), [[1.25, 3.25], [5.25, 7.25]])
        np.testing.assert_array_equal(dataset.GetRasterBand(indices["landcover"]).GetOverview(0).ReadAsArray(), [[10, 30], [60, 80]])
        np.testing.assert_array_equal(dataset.GetRasterBand(indices["woreda_id"]).GetOverview(0).ReadAsArray(), [[1, 3], [5, 7]])
        np.testing.assert_array_equal(dataset.GetRasterBand(indices["ethiopia_mask"]).GetOverview(0).ReadAsArray(), np.ones((2, 2)))

    def test_independent_reference_detects_average_used_for_categories(self):
        dataset = self.fixture()
        dataset.BuildOverviews("AVERAGE", [2])
        band = dataset.GetRasterBand(build.stack.band_indices_by_name(dataset)["landcover"])
        with mock.patch.object(build.grid, "COLS", 4), mock.patch.object(build.grid, "ROWS", 4):
            with self.assertRaisesRegex(ValueError, "method="):
                build.verify_overviews(band, band, {"name": "landcover", "overview_resample": "mode"}, [2], set(), {10, 20, 30, 40, 60, 80, 90, 95})

    def test_independent_mode_reference_passes(self):
        dataset = self.fixture()
        indices = build.stack.band_indices_by_name(dataset)
        build.apply_overview_policy(dataset, indices, [2])
        band = dataset.GetRasterBand(indices["landcover"])
        with mock.patch.object(build.grid, "COLS", 4), mock.patch.object(build.grid, "ROWS", 4):
            details = build.verify_overviews(band, band, {"name": "landcover", "overview_resample": "mode"}, [2], set(), {10, 20, 30, 40, 60, 80, 90, 95})
        self.assertEqual(details[0]["reference_mismatches"], 0)

    def test_lower_average_levels_preserve_nodata_contributor_weights(self):
        dataset = self.fixture()
        indices = build.stack.band_indices_by_name(dataset)
        values = np.full((4, 4), build.grid.NODATA, dtype=np.float32)
        values[:2, :2] = 1
        values[1, 1] = build.grid.NODATA
        values[0, 2] = 10
        values[2, 0] = 20
        values[2:, 2:] = 100
        band = dataset.GetRasterBand(indices["elevation_m"])
        band.WriteArray(values)
        build.apply_overview_policy(dataset, indices, [2, 4])
        actual = band.GetOverview(1).ReadAsArray()[0, 0]
        expected = np.float32(values[values != build.grid.NODATA].mean(dtype=np.float64))
        self.assertEqual(actual, expected)
        self.assertNotEqual(actual, np.float32((1 + 10 + 20 + 100) / 4))


class StatisticsTests(unittest.TestCase):
    def test_outliers_are_counted_not_modified(self):
        values = np.array([[5, 15, build.grid.NODATA]], dtype=np.float32)
        before = values.copy()
        result = build.band_statistics(values, np.ones_like(values, dtype=bool), {"name": "synthetic", "valid_min": 0, "valid_max": 10})
        self.assertEqual(result["above_sanity_range_cells"], 1)
        self.assertEqual(result["missing_inside_cells"], 1)
        np.testing.assert_array_equal(values, before)

    def test_unsupplied_bounds_are_null_not_invented(self):
        result = build.band_statistics(np.array([[5]], dtype=np.float32), np.array([[True]]), {"name": "synthetic", "valid_min": None, "valid_max": None})
        self.assertIsNone(result["below_sanity_range_cells"])
        self.assertIsNone(result["above_sanity_range_cells"])

    def test_nonfinite_output_fails(self):
        with self.assertRaisesRegex(ValueError, "non-finite"):
            build.band_statistics(np.array([[np.nan]], dtype=np.float32), np.array([[True]]), {"name": "synthetic", "valid_min": None, "valid_max": None})

    def test_json_missing_is_null_and_zero_is_not_missing(self):
        self.assertIsNone(build.json_value(build.grid.NODATA))
        self.assertIsNone(build.json_value(np.nan))
        self.assertEqual(build.json_value(0), 0.0)


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.prepared = self.root / "interim"
        self.processed = self.root / "processed"
        self.prepared.mkdir()
        self.processed.mkdir()
        self.paths = {name: self.prepared / filename for name, filename in (
            ("COG_PATH", "common_stack_v1.tif"), ("MANIFEST_PATH", "band_manifest.json"),
            ("LOOKUP_PATH", "woreda_id_lookup.csv"), ("VALIDATION_PATH", "stack_v1_validation.json"))}
        self.paths["COG_PATH"].write_bytes(b"synthetic raster bytes")
        self.manifest = {"build_script": "synthetic_build.py", "raster_sha256": build.sha256(self.paths["COG_PATH"])}
        self.paths["MANIFEST_PATH"].write_text(json.dumps(self.manifest))
        self.paths["LOOKUP_PATH"].write_text("synthetic lookup\n")
        self.paths["VALIDATION_PATH"].write_text('{"status": "passed"}\n')
        self.mocks = ExitStack()
        self.mocks.enter_context(mock.patch.object(build, "ROOT", self.root))
        self.mocks.enter_context(mock.patch.object(build, "PROCESSED", self.processed))
        for name, path in self.paths.items():
            self.mocks.enter_context(mock.patch.object(build, name, path))
        self.header = self.mocks.enter_context(mock.patch.object(build.stack, "assert_stack_matches_manifest"))

    def tearDown(self):
        self.mocks.close()
        self.temporary.cleanup()

    def test_finished_files_are_copied_identically(self):
        build.publish(self.manifest)
        for source in self.paths.values():
            self.assertEqual(source.read_bytes(), (self.processed / source.name).read_bytes())
        self.header.assert_called_once()

    def test_identical_owned_rerun_passes(self):
        build.publish(self.manifest)
        build.publish(self.manifest)
        self.assertEqual(self.header.call_count, 2)

    def test_unfamiliar_sidecar_is_not_overwritten_or_partially_published(self):
        destination = self.processed / self.paths["LOOKUP_PATH"].name
        destination.write_text("user-owned data\n")
        with self.assertRaisesRegex(ValueError, "unfamiliar deliverable"):
            build.publish(self.manifest)
        self.assertEqual(destination.read_text(), "user-owned data\n")
        self.assertFalse((self.processed / self.paths["COG_PATH"].name).exists())

    def test_existing_incomplete_pair_is_not_overwritten(self):
        destination = self.processed / self.paths["COG_PATH"].name
        destination.write_bytes(b"user-owned data")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            build.publish(self.manifest)
        self.assertEqual(destination.read_bytes(), b"user-owned data")

    def test_existing_different_version_content_is_not_overwritten(self):
        build.publish(self.manifest)
        destination = self.processed / self.paths["COG_PATH"].name
        destination.write_bytes(b"different raster")
        with self.assertRaisesRegex(ValueError, "explicit version decision"):
            build.publish(self.manifest)
        self.assertEqual(destination.read_bytes(), b"different raster")

    def test_unfamiliar_stack_builder_is_not_overwritten(self):
        build.publish(self.manifest)
        destination = self.processed / self.paths["MANIFEST_PATH"].name
        previous = copy.deepcopy(self.manifest)
        previous["build_script"] = "someone_else.py"
        destination.write_text(json.dumps(previous))
        with self.assertRaisesRegex(ValueError, "not built by this step"):
            build.publish(self.manifest)

    def test_legacy_builder_location_is_accepted_for_migration(self):
        (self.processed / self.paths["COG_PATH"].name).write_bytes(self.paths["COG_PATH"].read_bytes())
        destination = self.processed / self.paths["MANIFEST_PATH"].name
        previous = copy.deepcopy(self.manifest)
        previous["build_script"] = "data/interim/stack_v1/build_stack_v1.py"
        destination.write_text(json.dumps(previous))
        build.publish(self.manifest)
        self.assertEqual(json.loads(destination.read_text())["build_script"], self.manifest["build_script"])


if __name__ == "__main__":
    unittest.main()
