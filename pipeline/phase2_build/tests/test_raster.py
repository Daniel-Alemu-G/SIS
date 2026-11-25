import tempfile
import unittest
from pathlib import Path

import numpy as np
from osgeo import gdal

from config import stack, stack_v2
from pipeline.phase2_build import provenance, raster
from pipeline.phase2_build.tests.fixtures import calibrated_fixture, source_fixture
from scoring.engine import decode
from scoring.soil import calibrate_drought

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")


class RasterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="phase2-raster-", dir="/tmp/opencode")
        self.root = Path(self.temporary.name)
        self.source = source_fixture()
        self.indices = stack.band_indices_by_name(self.source)
        self.engine, self.manifest = calibrated_fixture(self.root, self.source)

    def tearDown(self):
        self.source = None
        self.temporary.cleanup()

    def test_windows_cover_each_pixel_once_including_edges(self):
        counts = np.zeros((9, 11), dtype=np.uint8)
        windows = list(raster.windows(11, 9, 4))
        self.assertEqual(windows[-1], (8, 8, 3, 1))
        for col, row, width, height in windows:
            counts[row:row + height, col:col + width] += 1
        np.testing.assert_array_equal(counts, np.ones_like(counts))

    def test_invalid_window_arguments_fail(self):
        for args in ((0, 9, 4), (11, -1, 4), (11, 9, True), (11, 9, 2.5)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                list(raster.windows(*args))

    def test_named_read_can_reorder_and_read_singletons(self):
        names = ("sg_phh2o", "ethiopia_mask")
        values = raster.read_named(self.source, self.indices, names, (1, 2, 1, 1))
        self.assertEqual(tuple(values), names)
        self.assertEqual(values["sg_phh2o"].shape, (1, 1))
        self.assertEqual(values["sg_phh2o"][0, 0], 6)
        singleton = raster.read_named(self.source, self.indices, ("sg_phh2o",))
        self.assertEqual(singleton["sg_phh2o"].shape, (9, 11))

    def test_missing_duplicate_and_empty_band_requests_fail(self):
        for names in ((), ("missing",), ("sg_phh2o", "sg_phh2o")):
            with self.subTest(names=names), self.assertRaises(ValueError):
                raster.read_named(self.source, self.indices, names)

    def test_non_float32_read_fails(self):
        source = gdal.GetDriverByName("MEM").Create("", 1, 1, 1, gdal.GDT_Byte)
        source.GetRasterBand(1).SetDescription("test")
        with self.assertRaisesRegex(ValueError, "float32"):
            raster.read_named(source, {"test": 1}, ("test",))

    def test_logical_missing_is_encoded_without_destroying_valid_zero(self):
        np.testing.assert_array_equal(raster.encode_score("best_crop", [0, 2, 255]), [0, 2, stack.grid.NODATA])
        np.testing.assert_array_equal(raster.encode_score("soil_quality", [0, np.nan, 1]), [0, stack.grid.NODATA, 1])
        with self.assertRaisesRegex(ValueError, "infinite"):
            raster.encode_score("soil_quality", [np.inf])

    def test_physical_contract_rejects_invalid_scores_and_classes(self):
        records = {row["name"]: row for row in self.manifest["bands"]}
        for name, values in (("soil_quality", [1.1]), ("best_crop", [255]), ("best_crop", [1.5]), ("ethiopia_mask", [stack.grid.NODATA])):
            with self.subTest(name=name, values=values), self.assertRaises(ValueError):
                raster.validate_physical(np.asarray(values, dtype=np.float32), records[name])
        with self.assertRaisesRegex(ValueError, "outside"):
            raster.validate_physical(np.asarray([0.1], dtype=np.float32), records["soil_quality"], inside=np.asarray([False]))

    def test_nonfinite_or_wrong_dtype_physical_values_fail(self):
        record = self.manifest["bands"][38]
        for values in (np.array([np.nan], dtype=np.float32), np.array([0.5], dtype=np.float64)):
            with self.assertRaisesRegex(ValueError, "finite float32"):
                raster.validate_physical(values, record)

    def test_window_calibration_equals_full_array(self):
        values = raster.read_named(self.source, self.indices, ("kiremt_mean_mm", "belg_mean_mm", "kiremt_cv", "belg_cv", "ethiopia_mask"))
        full = calibrate_drought(*(decode(values[name]) for name in ("kiremt_mean_mm", "belg_mean_mm", "kiremt_cv", "belg_cv")),
                                values["ethiopia_mask"] == 1, low_percentile=25, high_percentile=90, method="linear")
        for size in (1, 4, 127):
            self.assertEqual(raster.calibrate_stack(self.source, self.indices, self.engine.config.params, window_size=size), full)

    def test_bitwise_comparison_detects_negative_zero(self):
        self.assertEqual(raster.bitwise_mismatches(np.array([0], dtype=np.float32), np.array([-0.0], dtype=np.float32)), 1)
        with self.assertRaises(ValueError):
            raster.bitwise_mismatches(np.array([0]), np.array([0]))

    def test_cog_preserves_core_bits_tags_and_per_band_overviews(self):
        raw_path, cog_path = self.root / "raw.tif", self.root / "cog.tif"
        raster.write_stack(self.source, self.manifest["bands"], self.engine, raw_path, cog_path,
                           parent_sha256="1" * 64, overview_factors=[2, 4], window_size=3)
        cog = gdal.OpenEx(str(cog_path), gdal.OF_RASTER | gdal.OF_READONLY)
        raw = gdal.OpenEx(str(raw_path), gdal.OF_RASTER | gdal.OF_READONLY)
        indices = stack.band_indices_by_name(cog)
        self.assertEqual(len(indices), 60)
        self.assertEqual(cog.GetMetadata("IMAGE_STRUCTURE")["LAYOUT"], "COG")
        for name in stack.CORE_BAND_NAMES:
            original, copied = self.source.GetRasterBand(self.indices[name]), cog.GetRasterBand(indices[name])
            self.assertEqual(raster.bitwise_mismatches(original.ReadAsArray(), copied.ReadAsArray()), 0)
            self.assertEqual(original.GetMetadata(), copied.GetMetadata())
        checks = raster.verify_overviews(cog, self.manifest["bands"], [2, 4], raw=raw)
        self.assertEqual(len(checks), 60)
        self.assertTrue(all(level["reference_mismatches"] == 0 for check in checks for level in check["levels"]))
        self.assertEqual(next(check for check in checks if check["name"] == "best_crop")["levels"][0]["method"], "mode")
        cog = raw = None

    def test_different_window_sizes_produce_identical_cog_bytes(self):
        hashes = []
        for index, size in enumerate((3, 256)):
            path = self.root / f"cog-{index}.tif"
            raster.write_stack(self.source, self.manifest["bands"], self.engine, self.root / f"raw-{index}.tif", path,
                               parent_sha256="1" * 64, overview_factors=[2, 4], window_size=size)
            hashes.append(provenance.sha256(path))
        self.assertEqual(*hashes)

    def test_one_corrupted_overview_pixel_fails_independent_verification(self):
        raw_path, cog_path = self.root / "raw.tif", self.root / "cog.tif"
        raster.write_stack(self.source, self.manifest["bands"], self.engine, raw_path, cog_path,
                           parent_sha256="1" * 64, overview_factors=[2, 4])
        raw = gdal.OpenEx(str(raw_path), gdal.OF_RASTER | gdal.OF_UPDATE)
        indices = stack.band_indices_by_name(raw)
        raw.GetRasterBand(indices["soil_quality"]).GetOverview(1).WriteArray(np.asarray([[0.123]], dtype=np.float32), 1, 1)
        raw.FlushCache()
        with self.assertRaisesRegex(ValueError, "soil_quality overview 4.*mismatches"):
            raster.verify_overviews(raw, self.manifest["bands"], [2, 4])
        raw = None

    def test_assembly_cannot_overwrite_source_or_use_wrong_band_order(self):
        with self.assertRaisesRegex(ValueError, "in place"):
            raster.write_stack(self.source, self.manifest["bands"], self.engine, self.root / "same.tif", self.root / "same.tif",
                               parent_sha256="1" * 64, overview_factors=[])
        with self.assertRaisesRegex(ValueError, "named-band"):
            raster.write_stack(self.source, list(reversed(self.manifest["bands"])), self.engine, self.root / "raw.tif", self.root / "cog.tif",
                               parent_sha256="1" * 64, overview_factors=[])
