"""Synthetic checks for lossless, owner-approved phosphorus staging."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from osgeo import gdal, osr

import stage_phosphorus as stage


class PhosphorusTests(unittest.TestCase):
    def dataset(self, descriptions):
        dataset = gdal.GetDriverByName("MEM").Create(
            "", 1, 1, len(descriptions), gdal.GDT_Float32
        )
        for index, name in enumerate(descriptions, start=1):
            dataset.GetRasterBand(index).SetDescription(name)
        return dataset

    def test_owner_unit_is_preserved(self):
        self.assertEqual(stage.UNIT, "mg/100kg")
        values = np.array([[100, 1596.8274]], dtype=np.float32)
        staged, _ = stage.masked_copy(values, np.ones_like(values, dtype=bool))
        np.testing.assert_array_equal(staged.view(np.uint32), values.view(np.uint32))

    def test_source_gaps_and_nonfinite_values_stay_missing(self):
        values = np.array([[stage.NODATA, np.nan, np.inf, -np.inf, 5]], dtype=np.float32)
        staged, valid = stage.masked_copy(values, np.ones_like(values, dtype=bool))
        np.testing.assert_array_equal(valid, [[False, False, False, False, True]])
        np.testing.assert_array_equal(staged, [[stage.NODATA] * 4 + [5]])

    def test_outside_values_become_canonical_nodata(self):
        values = np.array([[100, 200]], dtype=np.float32)
        staged, _ = stage.masked_copy(values, np.array([[True, False]]))
        np.testing.assert_array_equal(staged, [[100, stage.NODATA]])

    def test_zero_is_not_reclassified_as_missing(self):
        staged, valid = stage.masked_copy(
            np.array([[0]], dtype=np.float32), np.array([[True]])
        )
        self.assertEqual(staged[0, 0], 0)
        self.assertTrue(valid[0, 0])

    def test_values_above_provisional_range_are_not_clipped(self):
        values = np.array([[9091.994]], dtype=np.float32)
        staged, _ = stage.masked_copy(values, np.array([[True]]))
        np.testing.assert_array_equal(staged, values)

    def test_mismatched_shapes_fail(self):
        with self.assertRaisesRegex(ValueError, "matching 2-D shapes"):
            stage.masked_copy(np.ones((2, 2)), np.ones((1, 2), dtype=bool))

    def test_nonboolean_mask_fails(self):
        with self.assertRaisesRegex(ValueError, "boolean country mask"):
            stage.masked_copy(np.ones((1, 1)), np.ones((1, 1)))

    def test_band_name_resolution_is_independent_of_order(self):
        dataset = self.dataset(["other", "sg_p"])
        band = stage.band_by_name(dataset, "sg_p")
        self.assertEqual(band.GetDescription(), "sg_p")

    def test_unnamed_singleton_requires_explicit_alias(self):
        dataset = self.dataset([""])
        with self.assertRaisesRegex(ValueError, "missing or duplicate"):
            stage.band_by_name(dataset, "sg_p")
        self.assertEqual(
            stage.band_by_name(dataset, "sg_p", unnamed_alias="sg_p").GetDescription(),
            "",
        )

    def test_alias_cannot_hide_an_unnamed_multiband_source(self):
        dataset = self.dataset(["", ""])
        with self.assertRaisesRegex(ValueError, "exactly one band"):
            stage.band_by_name(dataset, "sg_p", unnamed_alias="sg_p")

    def test_duplicate_names_fail(self):
        dataset = self.dataset(["sg_p", "sg_p"])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            stage.band_by_name(dataset, "sg_p")

    def test_unexpected_named_band_is_not_silently_aliased(self):
        dataset = self.dataset(["nitrogen"])
        with self.assertRaisesRegex(ValueError, "missing raster band"):
            stage.band_by_name(dataset, "sg_p", unnamed_alias="sg_p")

    def source_dataset(self, dtype=gdal.GDT_Float32):
        dataset = gdal.GetDriverByName("MEM").Create("", 1, 1, 1, dtype)
        band = dataset.GetRasterBand(1)
        band.SetNoDataValue(stage.NODATA)
        band.WriteArray(np.array([[100]], dtype=np.float32))
        return dataset

    def read_fixture(self, dataset):
        with (
            mock.patch.object(stage, "open_raster", return_value=dataset),
            mock.patch.object(stage, "COLS", 1),
            mock.patch.object(stage, "ROWS", 1),
            mock.patch.object(stage, "assert_matches_grid") as grid_check,
        ):
            result = stage.read_source()
            grid_check.assert_called_once_with(stage.SOURCE_PATH)
            return result

    def test_source_identity_and_owner_approval_are_reported_honestly(self):
        values, metadata = self.read_fixture(self.source_dataset())
        np.testing.assert_array_equal(values, [[100]])
        self.assertEqual(metadata["band_description"], "")
        self.assertEqual(metadata["upstream_product"], "AfricanSoilGrids")
        self.assertIn("none", metadata["required_scaling_and_back_transform"])
        self.assertEqual(metadata["provenance_authority"], stage.OWNER_APPROVAL)
        self.assertIsNone(metadata["product_version"])
        self.assertIsNone(metadata["embedded_unit"])

    def test_wrong_source_dtype_fails(self):
        with self.assertRaisesRegex(ValueError, "must be Float32"):
            self.read_fixture(self.source_dataset(gdal.GDT_Int16))

    def test_wrong_source_nodata_fails(self):
        dataset = self.source_dataset()
        dataset.GetRasterBand(1).SetNoDataValue(0)
        with self.assertRaisesRegex(ValueError, "nodata must equal canonical"):
            self.read_fixture(dataset)

    def test_unapplied_source_scale_fails(self):
        dataset = self.source_dataset()
        dataset.GetRasterBand(1).SetScale(0.01)
        with self.assertRaisesRegex(ValueError, "embedded scale/offset"):
            self.read_fixture(dataset)

    def test_unapplied_source_offset_fails(self):
        dataset = self.source_dataset()
        dataset.GetRasterBand(1).SetOffset(1)
        with self.assertRaisesRegex(ValueError, "embedded scale/offset"):
            self.read_fixture(dataset)

    def test_conflicting_embedded_unit_fails(self):
        dataset = self.source_dataset()
        dataset.GetRasterBand(1).SetUnitType("mg/kg")
        with self.assertRaisesRegex(ValueError, "disagrees"):
            self.read_fixture(dataset)

    def test_conflicting_unit_metadata_fails(self):
        dataset = self.source_dataset()
        dataset.GetRasterBand(1).SetMetadataItem("unit", "log-scaled")
        with self.assertRaisesRegex(ValueError, "disagrees"):
            self.read_fixture(dataset)

    def test_written_units_values_and_inclusion_approval_round_trip(self):
        with tempfile.TemporaryDirectory(prefix="p1-p-tests-", dir="/tmp/opencode") as directory:
            path = Path(directory) / "staged.tif"
            values = np.array([[100, stage.NODATA], [9091.994, 0]], dtype=np.float32)
            inside = np.array([[True, False], [True, True]])
            with (
                mock.patch.object(stage, "OUTPUT_TEMP_PATH", path),
                mock.patch.object(stage, "COLS", 2),
                mock.patch.object(stage, "ROWS", 2),
                mock.patch.object(stage, "assert_matches_grid") as grid_check,
            ):
                stage.write_output(values)
                checks = stage.validate_output(values, inside)
                grid_check.assert_called_once_with(path)
            self.assertEqual(checks["bitwise_round_trip_mismatches"], 0)
            self.assertEqual(checks["populated_outside_cells_after_mask"], 0)
            dataset = stage.open_raster(path)
            self.assertEqual(dataset.GetGeoTransform(), stage.TRANSFORM.to_gdal())
            expected_crs = osr.SpatialReference()
            expected_crs.SetFromUserInput(stage.CRS)
            actual_crs = dataset.GetSpatialRef().Clone()
            for spatial_ref in (actual_crs, expected_crs):
                spatial_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            self.assertTrue(actual_crs.IsSame(expected_crs))
            self.assertEqual(dataset.GetMetadataItem("AREA_OR_POINT"), "Area")
            self.assertEqual(dataset.GetMetadataItem("publication_ready"), "true")
            self.assertEqual(dataset.GetMetadataItem("stack_eligible"), "true")
            self.assertIn("AfricanSoilGrids", stage.band_by_name(dataset, "sg_p").GetMetadataItem("source"))
            self.assertEqual(stage.band_by_name(dataset, "sg_p").GetUnitType(), "mg/100kg")
            dataset = None


if __name__ == "__main__":
    unittest.main()
