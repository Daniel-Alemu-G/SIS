"""Small synthetic regression checks for the P1-10 staging step.

These fixtures are not geographical rasters and never read or modify project
inputs. The real grid and full source round trips are checked by the build.
"""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from osgeo import gdal

import stage_landcover_climate as stage


class StagingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="p1-10-tests-", dir="/tmp/opencode"
        )
        self.addCleanup(self.temporary.cleanup)
        self.here = Path(self.temporary.name)

    def dataset(self, descriptions, values=None):
        if values is None:
            values = [np.ones((1, 4), dtype=np.float32) for _ in descriptions]
        rows, cols = values[0].shape
        dataset = gdal.GetDriverByName("MEM").Create(
            "", cols, rows, len(descriptions), gdal.GDT_Float32
        )
        for index, (name, array) in enumerate(zip(descriptions, values), start=1):
            band = dataset.GetRasterBand(index)
            band.SetDescription(name)
            band.SetNoDataValue(stage.NODATA)
            band.WriteArray(array)
        return dataset

    def records(self):
        records = []
        for index, name in enumerate(stage.RAINFALL_NAMES, start=1):
            season, stat = name.split("_", 1)
            records.append(
                {
                    "band": index,
                    "name": name,
                    "season": season,
                    "stat": stat,
                    "months": stage.SEASON_MONTHS[season],
                    "years_used": list(range(2000, 2015)),
                    "n_years": 15,
                }
            )
        return copy.deepcopy(records)

    def load_records(self, records):
        path = self.here / "rainfall.json"
        path.write_text(json.dumps(records))
        with mock.patch.object(stage, "RAINFALL_METADATA_PATH", path):
            return stage.load_rainfall_metadata()

    def test_band_indices_resolve_names_after_reordering(self):
        dataset = self.dataset(["tmean_bega", "tmean_kiremt", "tmean_belg"])
        self.assertEqual(
            stage.band_indices_by_name(dataset),
            {"tmean_bega": 1, "tmean_kiremt": 2, "tmean_belg": 3},
        )

    def test_duplicate_names_fail(self):
        dataset = self.dataset(["rain", "rain"])
        with self.assertRaisesRegex(ValueError, "duplicate source band name"):
            stage.band_indices_by_name(dataset)

    def test_empty_names_fail(self):
        dataset = self.dataset([""])
        with self.assertRaisesRegex(ValueError, "missing or duplicate"):
            stage.band_indices_by_name(dataset)

    def test_valid_rainfall_metadata_passes(self):
        records = self.records()
        self.assertEqual(self.load_records(records), records)

    def test_rainfall_metadata_must_contain_nine_records(self):
        with self.assertRaisesRegex(ValueError, "nine records"):
            self.load_records(self.records()[:-1])

    def test_rainfall_metadata_record_must_be_an_object(self):
        records = self.records()
        records[0] = None
        with self.assertRaisesRegex(ValueError, "not an object"):
            self.load_records(records)

    def test_wrong_metadata_name_fails(self):
        records = self.records()
        records[0]["name"] = "belg_mean_mm"
        with self.assertRaisesRegex(ValueError, "metadata name"):
            self.load_records(records)

    def test_wrong_metadata_band_fails(self):
        records = self.records()
        records[0]["band"] = 2
        with self.assertRaisesRegex(ValueError, "metadata band"):
            self.load_records(records)

    def test_wrong_season_fails_even_with_matching_months(self):
        records = self.records()
        records[0]["season"] = "belg"
        records[0]["months"] = stage.SEASON_MONTHS["belg"]
        with self.assertRaisesRegex(ValueError, "approved"):
            self.load_records(records)

    def test_wrong_season_months_fail(self):
        records = self.records()
        records[0]["months"] = [6, 7, 8]
        with self.assertRaisesRegex(ValueError, "approved"):
            self.load_records(records)

    def test_duplicate_years_fail(self):
        records = self.records()
        records[0]["years_used"][1] = records[0]["years_used"][0]
        with self.assertRaisesRegex(ValueError, "year metadata"):
            self.load_records(records)

    def test_noninteger_years_fail(self):
        records = self.records()
        records[0]["years_used"][0] = "2000"
        with self.assertRaisesRegex(ValueError, "year metadata"):
            self.load_records(records)

    def test_unsorted_years_fail(self):
        records = self.records()
        records[0]["years_used"].reverse()
        with self.assertRaisesRegex(ValueError, "year metadata"):
            self.load_records(records)

    def test_stat_year_sets_must_match_within_a_season(self):
        records = self.records()
        records[1]["years_used"] = list(range(2001, 2016))
        with self.assertRaisesRegex(ValueError, "different season-years"):
            self.load_records(records)

    def test_wrong_stat_fails(self):
        records = self.records()
        records[0]["stat"] = "cv"
        with self.assertRaisesRegex(ValueError, "stat"):
            self.load_records(records)

    def test_temperature_mapping_uses_names_not_source_positions(self):
        datasets = {
            stage.LANDCOVER_PATH: self.dataset(["Map"]),
            stage.RAINFALL_PATH: self.dataset(stage.RAINFALL_NAMES),
            stage.TEMPERATURE_PATH: self.dataset(
                ["tmean_bega", "tmean_kiremt", "tmean_belg"]
            ),
        }
        with mock.patch.object(stage, "open_raster", side_effect=datasets.__getitem__):
            definitions = stage.build_band_definitions(self.records())
        temperatures = [d for d in definitions if d["group"] == "temperature"]
        self.assertEqual(
            [d["name"] for d in temperatures],
            ["tmean_kiremt_c", "tmean_bega_c", "tmean_belg_c"],
        )
        self.assertEqual([d["source_band"] for d in temperatures], [2, 1, 3])

    def rainfall_dataset(self, *, incorrect=False, all_missing=False):
        descriptions, values = [], []
        for season in stage.SEASON_MONTHS:
            mean = np.array([[10, 20, stage.NODATA, 40]], dtype=np.float32)
            standard_deviation = np.array(
                [[2, 4, stage.NODATA, 8]], dtype=np.float32
            )
            coefficient = np.array(
                [[0.2, 0.2, stage.NODATA, 0.2]], dtype=np.float32
            )
            if incorrect:
                coefficient[0, 0] = 0.3
            if all_missing:
                mean.fill(stage.NODATA)
                standard_deviation.fill(stage.NODATA)
                coefficient.fill(stage.NODATA)
            descriptions.extend(
                [f"{season}_mean_mm", f"{season}_std_mm", f"{season}_cv"]
            )
            values.extend([mean, standard_deviation, coefficient])
        return self.dataset(descriptions, values)

    def test_cv_consistency_passes_without_counting_nodata(self):
        dataset = self.rainfall_dataset()
        with mock.patch.object(stage, "open_raster", return_value=dataset):
            checks = stage.validate_rainfall_cv()
        self.assertTrue(all(c["compared_cells"] == 3 for c in checks.values()))
        self.assertTrue(all(c["tolerance_mismatches"] == 0 for c in checks.values()))

    def test_cv_consistency_rejects_incorrect_coefficients(self):
        dataset = self.rainfall_dataset(incorrect=True)
        with mock.patch.object(stage, "open_raster", return_value=dataset):
            with self.assertRaisesRegex(ValueError, "inconsistent with std/mean"):
                stage.validate_rainfall_cv()

    def test_cv_consistency_fails_clearly_when_no_values_can_be_compared(self):
        dataset = self.rainfall_dataset(all_missing=True)
        with mock.patch.object(stage, "open_raster", return_value=dataset):
            with self.assertRaisesRegex(ValueError, "no valid positive-mean cells"):
                stage.validate_rainfall_cv()

    def test_nodata_and_nonfinite_values_remain_missing(self):
        values = np.array([0, stage.NODATA, np.nan, np.inf, 10], dtype=np.float32)
        np.testing.assert_array_equal(
            stage.valid_mask(values, stage.NODATA), [True, False, False, False, True]
        )
        np.testing.assert_array_equal(
            stage.valid_mask(values, 0), [False, True, False, False, True]
        )


if __name__ == "__main__":
    unittest.main()
