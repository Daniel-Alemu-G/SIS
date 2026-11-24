#!/usr/bin/env python3
"""Derive canonical-grid mean slope degrees from native 90 m percentages.

The conversion is deliberately applied to every valid 90 m source pixel before
GDAL averages those degree values onto the common grid. The existing 1 km
percentage raster is checked for alignment but is never used as an input.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config.grid import (  # noqa: E402
    COLS,
    CRS,
    EAST,
    NODATA,
    NORTH,
    ROWS,
    SOUTH,
    WEST,
    assert_matches_grid,
)

SOURCE_PATH = ROOT / "data/inputs/Slope/slope_albers_90m.tif"
REFERENCE_PATH = ROOT / "data/inputs/Slope/slope_1km.tif"
MASK_PATH = ROOT / "data/interim/admin_grid/ethiopia_mask.tif"
HERE = ROOT / "data/interim/slope_mean_deg"
OUTPUT_PATH = HERE / "slope_mean_deg.tif"
REPORT_PATH = HERE / "slope_mean_deg_report.json"
DEGREES_TEMP_PATH = HERE / "slope_degrees_90m.tmp.tif"
WARP_TEMP_PATH = HERE / "slope_mean_deg.tmp.tif"
REPORT_TEMP_PATH = HERE / "slope_mean_deg_report.tmp.json"

BLOCK_ROWS = 256
VALID_MIN_DEGREES = 0.0
VALID_MAX_DEGREES = 60.0

gdal.UseExceptions()


def relative(path: Path) -> str:
    """Return a stable repository-relative path for reports."""
    return path.relative_to(ROOT).as_posix()


def open_raster(path: Path, access: int = gdal.GA_ReadOnly):
    """Open a raster or fail with its path in the diagnostic."""
    dataset = gdal.Open(str(path), access)
    if dataset is None:
        raise RuntimeError(f"cannot open raster: {relative(path)}")
    return dataset


def nodata_mask(array: np.ndarray, nodata: float | None) -> np.ndarray:
    """Return true where an array contains its declared nodata value."""
    if nodata is None:
        raise ValueError("source raster must declare nodata")
    if math.isnan(nodata):
        return np.isnan(array)
    return array == nodata


def convert_native_percent_to_degrees() -> dict:
    """Block-convert source percentages to degrees on the native 90 m grid."""
    source = open_raster(SOURCE_PATH)
    if source.RasterCount != 1:
        raise ValueError(
            f"{relative(SOURCE_PATH)} has {source.RasterCount} bands; expected 1"
        )
    if (source.GetMetadataItem("AREA_OR_POINT") or "").lower() != "area":
        raise ValueError(f"{relative(SOURCE_PATH)} is not tagged pixel-is-area")

    spatial_ref = source.GetSpatialRef()
    if spatial_ref is None or not spatial_ref.IsProjected():
        raise ValueError(
            f"{relative(SOURCE_PATH)} must retain its projected Albers CRS"
        )

    source_size = [source.RasterXSize, source.RasterYSize]
    source_transform = list(source.GetGeoTransform())
    source_crs_name = spatial_ref.GetName()
    source_projection_method = spatial_ref.GetAttrValue("PROJECTION")
    source_proj4 = spatial_ref.ExportToProj4()

    source_band = source.GetRasterBand(1)
    source_nodata = source_band.GetNoDataValue()
    if source_nodata is None:
        raise ValueError(f"{relative(SOURCE_PATH)} does not declare nodata")

    driver = gdal.GetDriverByName("GTiff")
    converted = driver.Create(
        str(DEGREES_TEMP_PATH),
        source.RasterXSize,
        source.RasterYSize,
        1,
        gdal.GDT_Float32,
        options=[
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
            "BIGTIFF=IF_SAFER",
            "NUM_THREADS=ALL_CPUS",
        ],
    )
    if converted is None:
        raise RuntimeError(f"cannot create {relative(DEGREES_TEMP_PATH)}")
    converted.SetGeoTransform(source.GetGeoTransform())
    converted.SetSpatialRef(spatial_ref)
    converted.SetMetadataItem("AREA_OR_POINT", "Area")
    converted.SetMetadataItem("native_unit", "percent")
    converted.SetMetadataItem("conversion", "degrees(arctan(percent / 100))")
    converted_band = converted.GetRasterBand(1)
    converted_band.SetNoDataValue(NODATA)
    converted_band.SetDescription("slope_deg_90m_intermediate")
    converted_band.SetMetadataItem("unit", "degrees")

    valid_count = 0
    tagged_nodata_count = 0
    nonfinite_count = 0
    source_min = math.inf
    source_max = -math.inf
    degree_min = math.inf
    degree_max = -math.inf
    degree_sum = 0.0

    for y_offset in range(0, source.RasterYSize, BLOCK_ROWS):
        height = min(BLOCK_ROWS, source.RasterYSize - y_offset)
        values = source_band.ReadAsArray(
            0, y_offset, source.RasterXSize, height
        )
        if values is None:
            raise RuntimeError(f"failed reading source rows at {y_offset}")

        finite = np.isfinite(values)
        tagged_nodata = nodata_mask(values, source_nodata)
        valid = finite & ~tagged_nodata
        negative_count = int(np.count_nonzero(valid & (values < 0)))
        if negative_count:
            raise ValueError(
                f"{relative(SOURCE_PATH)} has {negative_count} negative "
                f"non-nodata slope values in rows starting at {y_offset}"
            )

        output = np.full(values.shape, NODATA, dtype=np.float32)
        if np.any(valid):
            source_values = values[valid].astype(np.float64, copy=False)
            degree_values = np.degrees(np.arctan(source_values / 100.0))
            output[valid] = degree_values.astype(np.float32)

            count = int(source_values.size)
            valid_count += count
            source_min = min(source_min, float(source_values.min()))
            source_max = max(source_max, float(source_values.max()))
            degree_min = min(degree_min, float(degree_values.min()))
            degree_max = max(degree_max, float(degree_values.max()))
            degree_sum += float(degree_values.sum(dtype=np.float64))

        tagged_nodata_count += int(np.count_nonzero(tagged_nodata))
        nonfinite_count += int(np.count_nonzero(~finite & ~tagged_nodata))
        converted_band.WriteArray(output, 0, y_offset)

    if valid_count == 0:
        raise ValueError(f"{relative(SOURCE_PATH)} has no valid slope values")

    converted_band.FlushCache()
    converted.FlushCache()
    converted = None
    source = None

    return {
        "path": relative(SOURCE_PATH),
        "size": source_size,
        "crs": {
            "name": source_crs_name,
            "projection_method": source_projection_method,
            "proj4": source_proj4,
        },
        "geotransform": source_transform,
        "pixel_size_m": [abs(source_transform[1]), abs(source_transform[5])],
        "native_unit": "percent",
        "nodata": source_nodata,
        "valid_pixels": valid_count,
        "tagged_nodata_pixels": tagged_nodata_count,
        "nonfinite_pixels": nonfinite_count,
        "source_percent_min": source_min,
        "source_percent_max": source_max,
        "converted_degree_min": degree_min,
        "converted_degree_mean": degree_sum / valid_count,
        "converted_degree_max": degree_max,
    }


def average_degrees_to_common_grid() -> None:
    """Warp the converted 90 m values with average resampling."""
    options = gdal.WarpOptions(
        format="GTiff",
        dstSRS=CRS,
        outputBounds=(WEST, SOUTH, EAST, NORTH),
        width=COLS,
        height=ROWS,
        outputType=gdal.GDT_Float32,
        srcNodata=NODATA,
        dstNodata=NODATA,
        resampleAlg="average",
        multithread=True,
        warpOptions=["INIT_DEST=NO_DATA", "UNIFIED_SRC_NODATA=YES"],
        creationOptions=[
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
            "NUM_THREADS=ALL_CPUS",
        ],
    )
    output = gdal.Warp(
        str(WARP_TEMP_PATH), str(DEGREES_TEMP_PATH), options=options
    )
    if output is None:
        raise RuntimeError("GDAL Warp did not create the common-grid slope")
    output.SetMetadataItem("AREA_OR_POINT", "Area")
    output.SetMetadataItem("source", relative(SOURCE_PATH))
    output.SetMetadataItem("native_unit", "percent")
    output.SetMetadataItem("conversion", "degrees(arctan(percent / 100))")
    output.SetMetadataItem("resample", "average after per-pixel conversion")
    band = output.GetRasterBand(1)
    band.SetDescription("slope_mean_deg")
    band.SetMetadataItem("unit", "degrees")
    band.SetMetadataItem("source", relative(SOURCE_PATH))
    band.SetMetadataItem("resample", "average")
    output.FlushCache()
    output = None
    assert_matches_grid(WARP_TEMP_PATH)


def mask_and_summarize_output() -> dict:
    """Apply the approved Ethiopia mask and collect exact output statistics."""
    assert_matches_grid(MASK_PATH)
    mask = open_raster(MASK_PATH)
    output = open_raster(WARP_TEMP_PATH, gdal.GA_Update)
    mask_band = mask.GetRasterBand(1)
    output_band = output.GetRasterBand(1)

    inside_count = 0
    outside_count = 0
    valid_inside_count = 0
    missing_inside_count = 0
    populated_outside_before_mask = 0
    below_range_count = 0
    above_range_count = 0
    value_min = math.inf
    value_max = -math.inf
    value_sum = 0.0

    for y_offset in range(0, ROWS, BLOCK_ROWS):
        height = min(BLOCK_ROWS, ROWS - y_offset)
        mask_values = mask_band.ReadAsArray(0, y_offset, COLS, height)
        values = output_band.ReadAsArray(0, y_offset, COLS, height)
        if mask_values is None or values is None:
            raise RuntimeError(f"failed reading common-grid rows at {y_offset}")
        if not np.all(np.isin(mask_values, (0, 1))):
            raise ValueError("ethiopia_mask contains values other than 0 and 1")

        inside = mask_values == 1
        outside = ~inside
        finite = np.isfinite(values)
        valid = finite & (values != NODATA)
        valid_inside = inside & valid

        inside_count += int(np.count_nonzero(inside))
        outside_count += int(np.count_nonzero(outside))
        valid_inside_count += int(np.count_nonzero(valid_inside))
        missing_inside_count += int(np.count_nonzero(inside & ~valid))
        populated_outside_before_mask += int(np.count_nonzero(outside & valid))

        if np.any(valid_inside):
            inside_values = values[valid_inside].astype(np.float64, copy=False)
            value_min = min(value_min, float(inside_values.min()))
            value_max = max(value_max, float(inside_values.max()))
            value_sum += float(inside_values.sum(dtype=np.float64))
            below_range_count += int(
                np.count_nonzero(inside_values < VALID_MIN_DEGREES)
            )
            above_range_count += int(
                np.count_nonzero(inside_values > VALID_MAX_DEGREES)
            )

        values[~valid | outside] = NODATA
        output_band.WriteArray(values.astype(np.float32, copy=False), 0, y_offset)

    if valid_inside_count == 0:
        raise ValueError("derived slope has no valid cells inside Ethiopia")

    output_band.FlushCache()
    output.FlushCache()
    output = None
    mask = None

    assert_matches_grid(WARP_TEMP_PATH)
    return {
        "inside_mask_cells": inside_count,
        "outside_mask_cells": outside_count,
        "valid_inside_cells": valid_inside_count,
        "missing_inside_cells": missing_inside_count,
        "populated_outside_cells_before_mask": populated_outside_before_mask,
        "outside_cells_after_mask": "all nodata",
        "min_deg": value_min,
        "mean_deg": value_sum / valid_inside_count,
        "max_deg": value_max,
        "below_sanity_range_cells": below_range_count,
        "above_sanity_range_cells": above_range_count,
        "sanity_range_deg": [VALID_MIN_DEGREES, VALID_MAX_DEGREES],
    }


def compare_with_reference_percentage(derived_path: Path) -> dict:
    """Audit, but never derive from, the existing mean-percentage product."""
    assert_matches_grid(REFERENCE_PATH)
    assert_matches_grid(derived_path)
    reference = open_raster(REFERENCE_PATH)
    output = open_raster(derived_path)
    reference_band = reference.GetRasterBand(1)
    output_band = output.GetRasterBand(1)
    reference_nodata = reference_band.GetNoDataValue()

    compared_count = 0
    differing_count = 0
    absolute_difference_sum = 0.0
    maximum_absolute_difference = 0.0

    for y_offset in range(0, ROWS, BLOCK_ROWS):
        height = min(BLOCK_ROWS, ROWS - y_offset)
        percentages = reference_band.ReadAsArray(0, y_offset, COLS, height)
        correct_degrees = output_band.ReadAsArray(0, y_offset, COLS, height)
        if percentages is None or correct_degrees is None:
            raise RuntimeError(f"failed reading comparison rows at {y_offset}")

        reference_valid = np.isfinite(percentages) & ~nodata_mask(
            percentages, reference_nodata
        )
        output_valid = np.isfinite(correct_degrees) & (correct_degrees != NODATA)
        comparable = reference_valid & output_valid & (percentages >= 0)
        if not np.any(comparable):
            continue

        shortcut = np.degrees(
            np.arctan(percentages[comparable].astype(np.float64) / 100.0)
        )
        differences = np.abs(
            correct_degrees[comparable].astype(np.float64) - shortcut
        )
        compared_count += int(differences.size)
        differing_count += int(np.count_nonzero(differences > 1e-6))
        absolute_difference_sum += float(differences.sum(dtype=np.float64))
        maximum_absolute_difference = max(
            maximum_absolute_difference, float(differences.max())
        )

    reference = None
    output = None
    if compared_count == 0:
        raise ValueError("no cells can be compared with the 1 km reference")

    return {
        "path": relative(REFERENCE_PATH),
        "unit": "mean source percent",
        "grid_check": "passed",
        "role": "reference only; not used to derive slope_mean_deg",
        "incorrect_shortcut": "degrees(arctan(mean_percent / 100))",
        "compared_cells": compared_count,
        "cells_different_from_shortcut": differing_count,
        "mean_absolute_difference_deg": absolute_difference_sum / compared_count,
        "maximum_absolute_difference_deg": maximum_absolute_difference,
    }


def main() -> None:
    """Build the staged band and atomically publish its audit report."""
    for path in (DEGREES_TEMP_PATH, WARP_TEMP_PATH, REPORT_TEMP_PATH):
        path.unlink(missing_ok=True)

    try:
        source_summary = convert_native_percent_to_degrees()
        average_degrees_to_common_grid()
        output_summary = mask_and_summarize_output()
        reference_summary = compare_with_reference_percentage(WARP_TEMP_PATH)

        report = {
            "status": "passed",
            "output": relative(OUTPUT_PATH),
            "band": "slope_mean_deg",
            "unit": "degrees",
            "physical_dtype": "float32",
            "nodata": NODATA,
            "grid_check": "passed",
            "derivation": {
                "formula": "degrees(arctan(percent / 100))",
                "order": "convert each valid 90 m pixel, then average",
                "resampling": "GDAL average from projected 90 m grid",
                "masking": "set every cell outside ethiopia_mask to nodata",
            },
            "source_90m": source_summary,
            "output_statistics": output_summary,
            "existing_1km_percentage": reference_summary,
            "sources_unchanged": True,
        }
        REPORT_TEMP_PATH.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        WARP_TEMP_PATH.replace(OUTPUT_PATH)
        assert_matches_grid(OUTPUT_PATH)
        REPORT_TEMP_PATH.replace(REPORT_PATH)
        print(json.dumps(report, indent=2, sort_keys=True))
    finally:
        DEGREES_TEMP_PATH.unlink(missing_ok=True)
        WARP_TEMP_PATH.unlink(missing_ok=True)
        REPORT_TEMP_PATH.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
