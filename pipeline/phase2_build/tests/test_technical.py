import tempfile
import unittest
from pathlib import Path

import numpy as np
from osgeo import gdal

from config import stack
from pipeline.phase2_build import raster, technical
from pipeline.phase2_build.tests.fixtures import calibrated_fixture, source_fixture


class TechnicalTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="phase2-technical-", dir="/tmp/opencode")
        self.root = Path(self.temporary.name)
        self.source = source_fixture()
        self.engine, self.manifest = calibrated_fixture(self.root, self.source)
        raw_path, cog_path = self.root / "raw.tif", self.root / "cog.tif"
        raster.write_stack(self.source, self.manifest["bands"], self.engine, raw_path, cog_path,
                           parent_sha256="1" * 64, overview_factors=[2, 4], window_size=4)
        stored = gdal.OpenEx(str(cog_path), gdal.OF_RASTER | gdal.OF_READONLY)
        self.candidate = gdal.Translate("", stored, format="MEM")
        stored = None
        self.indices = stack.band_indices_by_name(self.candidate)

    def tearDown(self):
        self.source = self.candidate = None
        self.temporary.cleanup()

    def change(self, name, value):
        self.candidate.GetRasterBand(self.indices[name]).WriteArray(np.asarray([[value]], dtype=np.float32), 3, 3)

    def test_all_pixels_match_v1_and_shared_engine_with_different_windows(self):
        bands, samples = technical.compare_pixels(self.source, self.candidate, self.manifest, window_size=3)
        self.assertEqual(len(bands), 60)
        self.assertEqual(samples["inside_cells"], 80)
        self.assertEqual(samples["outside_cells"], 19)
        for band in bands:
            self.assertEqual(band["valid_inside_cells"] + band["missing_inside_cells"], 80)
        self.assertGreater(samples["matched_cells"]["hard_zero_slope"], 0)
        self.assertGreater(samples["matched_cells"]["missing_basic_input"], 0)

    def test_at_least_20_representative_and_random_scalar_reports_match(self):
        _, samples = technical.compare_pixels(self.source, self.candidate, self.manifest)
        result = technical.scalar_checks(self.source, self.candidate, self.manifest, samples["representatives"])
        self.assertGreaterEqual(result["compared_cells"], 32)
        self.assertEqual(result["bitwise_mismatches"], 0)
        self.assertEqual(result, technical.scalar_checks(self.source, self.candidate, self.manifest, samples["representatives"]))
        excluded = next(cell for cell in result["cells"] if "built_up" in cell["selection"])
        self.assertIsNone(excluded["scores"]["viability_teff"])
        self.assertIsNone(excluded["scores"]["best_crop"])

    def test_any_changed_v1_pixel_fails(self):
        self.change("sg_p", 123)
        with self.assertRaisesRegex(ValueError, "sg_p.*pixel mismatches"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_any_changed_score_pixel_fails(self):
        self.change("viability_teff", 0.123)
        with self.assertRaisesRegex(ValueError, "viability_teff.*pixel mismatches"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_categorical_255_is_never_a_valid_physical_id(self):
        self.change("best_crop", 255)
        with self.assertRaisesRegex(ValueError, "illegal categorical"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_changed_v1_provenance_tag_fails(self):
        self.candidate.GetRasterBand(self.indices["sg_p"]).SetMetadataItem("opaque_original_provenance", "changed")
        with self.assertRaisesRegex(ValueError, "per-band metadata changed"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_changed_calibration_fails_before_scoring(self):
        self.manifest["scoring"]["calibration"]["cv_high"] += 0.01
        with self.assertRaisesRegex(ValueError, "fresh input-window calibration"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_outside_non_nodata_and_nonfinite_physical_values_fail(self):
        self.candidate.GetRasterBand(self.indices["soil_quality"]).WriteArray(np.asarray([[0.1]], dtype=np.float32), 0, 0)
        with self.assertRaisesRegex(ValueError, "outside Ethiopia"):
            technical.compare_pixels(self.source, self.candidate, self.manifest)

    def test_raw_assembly_pixel_mismatch_fails(self):
        raw = gdal.Translate("", self.candidate, format="MEM")
        raw.GetRasterBand(self.indices["sg_p"]).WriteArray(np.asarray([[123]], dtype=np.float32), 3, 3)
        with self.assertRaisesRegex(ValueError, "raw assembly"):
            technical.compare_pixels(self.source, self.candidate, self.manifest, raw=raw)
        raw = None
