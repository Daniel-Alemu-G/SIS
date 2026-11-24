#!/usr/bin/env python3
"""Prove the approved float32/-9999 stack and per-band overview policy.

This is deliberately synthetic. It writes only to its data/interim folder and
proves that a COG translation retains average overviews for continuous values
and mode overviews for categorical, identifier, and mask values.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config.grid import CRS, NODATA  # noqa: E402

HERE = ROOT / "data/interim/storage_prototype"
SOURCE_PATH = HERE / "prototype_source.tif"
OUTPUT_PATH = HERE / "prototype_stack_cog.tif"
REPORT_PATH = HERE / "prototype_report.json"
TEMP_PATH = HERE / "prototype_stack_cog.tmp.tif"

SIZE = 512
OVERVIEW_FACTORS = [2, 4]
BANDS = (
    ("continuous", "float32", "AVERAGE"),
    ("category", "uint8", "MODE"),
    ("identifier", "uint16", "MODE"),
    ("mask", "uint8", "MODE"),
)

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")


def remove_own_output(path: Path) -> None:
    if path.exists():
        path.unlink()


def synthetic_arrays() -> list[np.ndarray]:
    rows, cols = np.indices((SIZE, SIZE))

    continuous = (rows * 0.25 + cols * 0.5).astype("float32")
    category = np.ones((SIZE, SIZE), dtype="float32")
    category[(rows % 2 == 1) & (cols % 2 == 1)] = 2
    identifier = np.full((SIZE, SIZE), 100, dtype="float32")
    identifier[(rows % 2 == 1) & (cols % 2 == 1)] = 200
    mask = np.ones((SIZE, SIZE), dtype="float32")
    mask[(rows % 2 == 1) & (cols % 2 == 1)] = 0

    # Exercise the shared physical sentinel without putting nodata in the mask.
    for array in (continuous, category, identifier):
        array[:32, :32] = NODATA
    return [continuous, category, identifier, mask]


def create_source() -> None:
    remove_own_output(SOURCE_PATH)
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(
        str(SOURCE_PATH),
        SIZE,
        SIZE,
        len(BANDS),
        gdal.GDT_Float32,
        options=(
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
        ),
    )
    dataset.SetGeoTransform((32.9, 1 / 120, 0, 14.9, 0, -1 / 120))
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(CRS)
    dataset.SetSpatialRef(spatial_ref)
    dataset.SetMetadataItem("AREA_OR_POINT", "Area")
    dataset.SetMetadataItem("stack_version", "prototype")

    for index, (array, (name, logical_dtype, resampling)) in enumerate(
        zip(synthetic_arrays(), BANDS, strict=True), start=1
    ):
        band = dataset.GetRasterBand(index)
        band.WriteArray(array)
        band.SetNoDataValue(NODATA)
        band.SetDescription(name)
        band.SetMetadataItem("logical_dtype", logical_dtype)
        band.SetMetadataItem("overview_resampling", resampling.lower())

    dataset.FlushCache()
    if dataset.BuildOverviews("AVERAGE", OVERVIEW_FACTORS) != gdal.CE_None:
        raise RuntimeError("could not create prototype overview levels")

    for band_index in range(2, len(BANDS) + 1):
        band = dataset.GetRasterBand(band_index)
        for overview_index in range(band.GetOverviewCount()):
            if (
                gdal.RegenerateOverview(
                    band, band.GetOverview(overview_index), "MODE"
                )
                != gdal.CE_None
            ):
                raise RuntimeError(
                    f"could not regenerate mode overview for band {band_index}"
                )
    dataset.FlushCache()
    dataset = None


def create_cog() -> None:
    remove_own_output(TEMP_PATH)
    remove_own_output(OUTPUT_PATH)
    result = gdal.Translate(
        str(TEMP_PATH),
        str(SOURCE_PATH),
        format="COG",
        creationOptions=(
            "COMPRESS=DEFLATE",
            "PREDICTOR=FLOATING_POINT",
            "OVERVIEW_COMPRESS=DEFLATE",
            "OVERVIEW_PREDICTOR=FLOATING_POINT",
            "BLOCKSIZE=256",
            "OVERVIEWS=FORCE_USE_EXISTING",
            "NUM_THREADS=ALL_CPUS",
        ),
    )
    if result is None:
        raise RuntimeError("COG translation failed")
    result = None
    TEMP_PATH.replace(OUTPUT_PATH)


def validate() -> dict:
    source = gdal.OpenEx(str(SOURCE_PATH), gdal.OF_RASTER | gdal.OF_READONLY)
    output = gdal.OpenEx(str(OUTPUT_PATH), gdal.OF_RASTER | gdal.OF_READONLY)
    if source is None or output is None:
        raise RuntimeError("could not reopen prototype outputs")

    checks = {
        "layout_is_cog": output.GetMetadataItem("LAYOUT", "IMAGE_STRUCTURE")
        == "COG",
        "physical_type_is_float32": all(
            output.GetRasterBand(index).DataType == gdal.GDT_Float32
            for index in range(1, output.RasterCount + 1)
        ),
        "nodata_is_minus_9999": all(
            output.GetRasterBand(index).GetNoDataValue() == NODATA
            for index in range(1, output.RasterCount + 1)
        ),
        "all_overviews_preserved": True,
        "overview_pixels_preserved": True,
        "categorical_overviews_are_legal": True,
        "mask_has_no_nodata_pixels": True,
    }
    overview_shapes = {}
    overview_unique_values = {}
    legal_values = {2: {1.0, 2.0, float(NODATA)}, 3: {100.0, 200.0, float(NODATA)}, 4: {0.0, 1.0}}

    for band_index, (name, _, _) in enumerate(BANDS, start=1):
        source_band = source.GetRasterBand(band_index)
        output_band = output.GetRasterBand(band_index)
        checks["all_overviews_preserved"] &= (
            source_band.GetOverviewCount()
            == output_band.GetOverviewCount()
            == len(OVERVIEW_FACTORS)
        )
        overview_shapes[name] = []
        overview_unique_values[name] = []
        for overview_index in range(output_band.GetOverviewCount()):
            source_array = source_band.GetOverview(overview_index).ReadAsArray()
            output_array = output_band.GetOverview(overview_index).ReadAsArray()
            checks["overview_pixels_preserved"] &= np.array_equal(
                source_array, output_array
            )
            overview_shapes[name].append(list(output_array.shape))
            if band_index > 1:
                unique = {float(value) for value in np.unique(output_array)}
                overview_unique_values[name].append(sorted(unique))
                checks["categorical_overviews_are_legal"] &= unique <= legal_values[
                    band_index
                ]

    mask = output.GetRasterBand(4).ReadAsArray()
    checks["mask_has_no_nodata_pixels"] = not np.any(mask == NODATA)
    if not all(checks.values()):
        raise RuntimeError(f"storage prototype failed: {checks}")

    result = {
        "status": "passed",
        "policy": {
            "physical_dtype": "float32",
            "nodata": NODATA,
            "continuous_overviews": "average",
            "categorical_id_mask_overviews": "mode",
        },
        "checks": checks,
        "overview_shapes": overview_shapes,
        "overview_unique_values": overview_unique_values,
        "source": str(SOURCE_PATH.relative_to(ROOT)),
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
    }
    source = None
    output = None
    return result


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    create_source()
    create_cog()
    report = validate()
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
