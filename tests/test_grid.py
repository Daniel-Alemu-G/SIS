import unittest

from osgeo import gdal, osr

from config.grid import (
    COLS,
    CRS,
    ROWS,
    TRANSFORM,
    GridMismatchError,
    assert_matches_grid,
    grid_mismatches,
)

gdal.UseExceptions()


class GridValidationTests(unittest.TestCase):
    def make_dataset(
        self,
        *,
        cols=COLS,
        rows=ROWS,
        transform=None,
        crs=CRS,
        area_or_point="Area",
    ):
        dataset = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
        dataset.SetGeoTransform(transform or TRANSFORM.to_gdal())
        spatial_ref = osr.SpatialReference()
        spatial_ref.SetFromUserInput(crs)
        dataset.SetSpatialRef(spatial_ref)
        if area_or_point is not None:
            dataset.SetMetadataItem("AREA_OR_POINT", area_or_point)
        return dataset

    def test_matching_dataset_passes(self):
        dataset = self.make_dataset()

        self.assertEqual(grid_mismatches(dataset), [])
        self.assertIsNone(assert_matches_grid(dataset))

    def test_shape_mismatch_fails(self):
        dataset = self.make_dataset(cols=COLS + 1)

        with self.assertRaisesRegex(GridMismatchError, r"1813x1380"):
            assert_matches_grid(dataset)

    def test_rotation_mismatch_fails(self):
        transform = list(TRANSFORM.to_gdal())
        transform[2] = 0.001
        dataset = self.make_dataset(transform=tuple(transform))

        with self.assertRaisesRegex(GridMismatchError, r"transform\[2\]"):
            assert_matches_grid(dataset)

    def test_non_equivalent_crs_fails(self):
        dataset = self.make_dataset(crs="EPSG:3857")

        with self.assertRaisesRegex(GridMismatchError, r"crs:"):
            assert_matches_grid(dataset)

    def test_pixel_is_point_fails(self):
        dataset = self.make_dataset(area_or_point="Point")

        with self.assertRaisesRegex(GridMismatchError, r"AREA_OR_POINT='Point'"):
            assert_matches_grid(dataset)

    def test_missing_pixel_convention_fails(self):
        dataset = self.make_dataset(area_or_point=None)

        with self.assertRaisesRegex(GridMismatchError, r"AREA_OR_POINT=None"):
            assert_matches_grid(dataset)


if __name__ == "__main__":
    unittest.main()
