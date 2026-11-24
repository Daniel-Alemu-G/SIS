import unittest
from unittest import mock

from osgeo import gdal

from config import stack

gdal.UseExceptions()


def manifest_fixture():
    return {
        "stack_version": stack.STACK_VERSION,
        "grid": stack.grid_metadata(),
        "overview_factors": stack.overview_factors(stack.grid.COLS, stack.grid.ROWS),
        "bands": [{
            "index": index, "name": name, "description": f"Synthetic {name}",
            "unit": stack.CORE_UNITS[name], "physical_dtype": "float32",
            "logical_dtype": stack.logical_dtype(name), "nodata": stack.grid.NODATA,
            "overview_resample": stack.overview_method(name),
            "source": "AfricanSoilGrids" if name == "sg_p" else "Synthetic fixture",
            "source_file": "synthetic.tif", "resample": "none",
            "missing_value_semantics": "Synthetic fixture only",
            "transform_applied": "none",
        } for index, name in enumerate(stack.CORE_BAND_NAMES, start=1)],
    }


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.manifest = manifest_fixture()

    def test_core_has_38_unique_names(self):
        self.assertEqual(len(stack.CORE_BAND_NAMES), 38)
        self.assertEqual(len(set(stack.CORE_BAND_NAMES)), 38)

    def test_valid_manifest_passes(self):
        self.assertEqual(tuple(stack.validate_manifest(self.manifest)), stack.CORE_BAND_NAMES)

    def test_version_mismatch_fails(self):
        self.manifest["stack_version"] = 2
        with self.assertRaisesRegex(stack.StackContractError, "stack_version"):
            stack.validate_manifest(self.manifest)

    def test_boolean_version_fails(self):
        self.manifest["stack_version"] = True
        with self.assertRaisesRegex(stack.StackContractError, "stack_version"):
            stack.validate_manifest(self.manifest)

    def test_grid_mismatch_fails(self):
        self.manifest["grid"]["cols"] += 1
        with self.assertRaisesRegex(stack.StackContractError, "manifest grid"):
            stack.validate_manifest(self.manifest)

    def test_missing_band_fails(self):
        self.manifest["bands"].pop()
        with self.assertRaisesRegex(stack.StackContractError, "38 core bands"):
            stack.validate_manifest(self.manifest)

    def test_reordered_band_fails(self):
        self.manifest["bands"][0], self.manifest["bands"][1] = self.manifest["bands"][1], self.manifest["bands"][0]
        with self.assertRaisesRegex(stack.StackContractError, "must be 'ethiopia_mask'"):
            stack.validate_manifest(self.manifest)

    def test_wrong_index_fails(self):
        self.manifest["bands"][0]["index"] = 2
        with self.assertRaisesRegex(stack.StackContractError, "invalid manifest index"):
            stack.validate_manifest(self.manifest)

    def test_mask_keeps_logical_uint8_and_mode(self):
        mask = stack.validate_manifest(self.manifest)["ethiopia_mask"]
        self.assertEqual(mask["logical_dtype"], "uint8")
        self.assertEqual(mask["physical_dtype"], "float32")
        self.assertEqual(mask["overview_resample"], "mode")

    def test_categorical_average_overviews_fail(self):
        self.manifest["bands"][0]["overview_resample"] = "average"
        with self.assertRaisesRegex(stack.StackContractError, "overview_resample"):
            stack.validate_manifest(self.manifest)

    def test_phosphorus_cannot_be_rescaled_or_relabelled(self):
        phosphorus = next(b for b in self.manifest["bands"] if b["name"] == "sg_p")
        phosphorus["unit"] = "mg/kg"
        with self.assertRaisesRegex(stack.StackContractError, "unit must be 'mg/100kg'"):
            stack.validate_manifest(self.manifest)
        phosphorus["unit"] = "mg/100kg"
        phosphorus["transform_applied"] = "divide by 100"
        with self.assertRaisesRegex(stack.StackContractError, "without scaling"):
            stack.validate_manifest(self.manifest)

    def test_phosphorus_source_is_not_global_soilgrids(self):
        phosphorus = next(b for b in self.manifest["bands"] if b["name"] == "sg_p")
        phosphorus["source"] = "SoilGrids 2.0 (ISRIC)"
        with self.assertRaisesRegex(stack.StackContractError, "AfricanSoilGrids"):
            stack.validate_manifest(self.manifest)

    def test_unresolved_unit_or_source_fails(self):
        self.manifest["bands"][0]["source"] = "TBD"
        with self.assertRaisesRegex(stack.StackContractError, "unresolved source"):
            stack.validate_manifest(self.manifest)

    def test_missing_physical_nodata_fails(self):
        self.manifest["bands"][0]["nodata"] = None
        with self.assertRaisesRegex(stack.StackContractError, "nodata"):
            stack.validate_manifest(self.manifest)


class RasterHeaderTests(unittest.TestCase):
    def setUp(self):
        self.manifest = manifest_fixture()
        self.dataset = gdal.GetDriverByName("MEM").Create("", 1, 1, 38, gdal.GDT_Float32)
        self.dataset.SetMetadataItem("stack_version", str(stack.STACK_VERSION))
        for record in self.manifest["bands"]:
            band = self.dataset.GetRasterBand(record["index"])
            band.SetDescription(record["name"])
            band.SetNoDataValue(stack.grid.NODATA)
            band.SetUnitType(record["unit"])
            for key in ("unit", "source", "source_file", "logical_dtype", "overview_resample"):
                band.SetMetadataItem(key, record[key])

    def validate(self, *, require_cog=False):
        # These tiny fixtures isolate tags/version validation. Grid validation
        # is tested against real grid-sized datasets in test_grid.py and the
        # complete built raster in the acceptance script.
        with mock.patch.object(stack, "assert_matches_grid") as grid_check:
            result = stack.assert_stack_matches_manifest(self.dataset, self.manifest, require_cog=require_cog)
            grid_check.assert_called_once_with(self.dataset)
            return result

    def test_matching_header_passes(self):
        self.assertEqual(tuple(self.validate()), stack.CORE_BAND_NAMES)

    def test_raster_version_mismatch_fails(self):
        self.dataset.SetMetadataItem("stack_version", "2")
        with self.assertRaisesRegex(stack.StackContractError, "does not match manifest"):
            self.validate()

    def test_duplicate_raster_name_fails(self):
        self.dataset.GetRasterBand(2).SetDescription("ethiopia_mask")
        with self.assertRaisesRegex(stack.StackContractError, "duplicate"):
            self.validate()

    def test_raster_unit_mismatch_fails(self):
        self.dataset.GetRasterBand(1).SetUnitType("wrong")
        with self.assertRaisesRegex(stack.StackContractError, "raster unit"):
            self.validate()

    def test_raster_source_mismatch_fails(self):
        self.dataset.GetRasterBand(1).SetMetadataItem("source", "wrong")
        with self.assertRaisesRegex(stack.StackContractError, "raster source"):
            self.validate()

    def test_plain_raster_is_not_accepted_as_cog(self):
        with self.assertRaisesRegex(stack.StackContractError, "predictor-3 COG"):
            self.validate(require_cog=True)

    def test_name_lookup_handles_reordering(self):
        dataset = gdal.GetDriverByName("MEM").Create("", 1, 1, 2, gdal.GDT_Float32)
        dataset.GetRasterBand(1).SetDescription("sg_p")
        dataset.GetRasterBand(2).SetDescription("sg_n")
        self.assertEqual(stack.band_indices_by_name(dataset), {"sg_p": 1, "sg_n": 2})


if __name__ == "__main__":
    unittest.main()
