#!/usr/bin/env python3
"""Apply the approved EthioSIS grid/storage normalization without resampling.

The source has one empty column east of the canonical grid. This script first
verifies that column, then window-copies the canonical 1812 columns and maps
non-finite source values to the approved physical -9999 sentinel. Individual
nutrients remain independently missing even when confidence is positive.
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
    NODATA,
    ROWS,
    assert_matches_grid,
)

SOURCE_PATH = ROOT / "data/inputs/EthioSIS/ethiosis_IDW_estimates.tif"
HERE = ROOT / "data/interim/input_normalization"
OUTPUT_PATH = HERE / "ethiosis_IDW_estimates_grid_v1.tif"
TEMP_PATH = HERE / "ethiosis_IDW_estimates_grid_v1.tmp.tif"
REPORT_PATH = HERE / "normalization_report.json"

EXPECTED_BANDS = (
    "ca",
    "cu",
    "fe",
    "mg",
    "mn",
    "mo",
    "k",
    "si",
    "na",
    "s",
    "zn",
    "confidence",
)
CHUNK_ROWS = 256

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")


def valid_mask(array: np.ndarray, nodata: float | None) -> np.ndarray:
    valid = np.isfinite(array)
    if nodata is not None and math.isfinite(nodata):
        valid &= array != nodata
    return valid


def remove_own_output(path: Path) -> None:
    if path.exists():
        path.unlink()


def verify_source(source) -> tuple[dict, np.ndarray]:
    if source.RasterXSize != COLS + 1 or source.RasterYSize != ROWS:
        raise RuntimeError(
            "expected the inventoried 1813x1380 source before applying the "
            f"approved window; got {source.RasterXSize}x{source.RasterYSize}"
        )
    descriptions = tuple(
        source.GetRasterBand(index).GetDescription()
        for index in range(1, source.RasterCount + 1)
    )
    if descriptions != EXPECTED_BANDS:
        raise RuntimeError(
            f"unexpected EthioSIS band descriptions: {descriptions!r}"
        )

    # A VRT view proves that CRS, transform, dimensions and pixel convention
    # match after the approved window, without creating or resampling data.
    canonical_view = gdal.Translate(
        "", source, format="VRT", srcWin=(0, 0, COLS, ROWS)
    )
    assert_matches_grid(canonical_view)
    canonical_view = None

    discarded = {}
    for index, name in enumerate(EXPECTED_BANDS, start=1):
        band = source.GetRasterBand(index)
        array = band.ReadAsArray(COLS, 0, 1, ROWS)
        valid = valid_mask(array, band.GetNoDataValue())
        finite = np.isfinite(array)
        discarded[name] = {
            "finite_cells": int(finite.sum()),
            "valid_cells": int(valid.sum()),
            "nonzero_cells": int(np.count_nonzero(np.nan_to_num(array, nan=0))),
        }
        if name == "confidence":
            if not np.all(finite) or np.any(array != 0):
                raise RuntimeError(
                    "discarded confidence column is not entirely finite zero"
                )
        elif np.any(valid):
            raise RuntimeError(
                f"discarded eastern column contains valid {name} values"
            )

    confidence_band = source.GetRasterBand(EXPECTED_BANDS.index("confidence") + 1)
    confidence = confidence_band.ReadAsArray(0, 0, COLS, ROWS)
    confidence_valid = valid_mask(confidence, confidence_band.GetNoDataValue())
    if not np.all(confidence_valid):
        raise RuntimeError("canonical confidence window contains missing values")
    if np.any((confidence < 0) | (confidence > 1)):
        raise RuntimeError("confidence values fall outside 0..1")
    return discarded, confidence


def create_output(source, confidence: np.ndarray) -> list[dict]:
    remove_own_output(TEMP_PATH)
    driver = gdal.GetDriverByName("GTiff")
    output = driver.Create(
        str(TEMP_PATH),
        COLS,
        ROWS,
        source.RasterCount,
        gdal.GDT_Float32,
        options=(
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
            "BIGTIFF=IF_SAFER",
            "NUM_THREADS=ALL_CPUS",
        ),
    )
    output.SetGeoTransform(source.GetGeoTransform())
    output.SetSpatialRef(source.GetSpatialRef())
    output.SetMetadata(source.GetMetadata())
    output.SetMetadataItem("AREA_OR_POINT", "Area")
    output.SetMetadataItem("normalization", "approved_first_1812_columns")
    output.SetMetadataItem("physical_dtype", "float32")
    output.SetMetadataItem("physical_nodata", str(NODATA))

    band_reports = []
    for index, name in enumerate(EXPECTED_BANDS, start=1):
        source_band = source.GetRasterBand(index)
        output_band = output.GetRasterBand(index)
        output_band.SetDescription(name)
        output_band.SetNoDataValue(NODATA)
        output_band.SetMetadata(source_band.GetMetadata())
        output_band.SetMetadataItem("logical_dtype", "float32")
        output_band.SetMetadataItem("source_band", name)

        valid_count = missing_count = missing_with_positive_confidence = 0
        valid_with_zero_confidence = 0
        for y_offset in range(0, ROWS, CHUNK_ROWS):
            height = min(CHUNK_ROWS, ROWS - y_offset)
            source_array = source_band.ReadAsArray(0, y_offset, COLS, height)
            valid = valid_mask(source_array, source_band.GetNoDataValue())
            normalized = np.where(valid, source_array, NODATA).astype(
                "float32", copy=False
            )
            output_band.WriteArray(normalized, 0, y_offset)

            valid_count += int(valid.sum())
            missing_count += int(valid.size - valid.sum())
            if name != "confidence":
                confidence_chunk = confidence[y_offset : y_offset + height]
                missing_with_positive_confidence += int(
                    ((confidence_chunk > 0) & ~valid).sum()
                )
                valid_with_zero_confidence += int(
                    ((confidence_chunk == 0) & valid).sum()
                )

        if name != "confidence" and valid_with_zero_confidence:
            raise RuntimeError(
                f"{name} has {valid_with_zero_confidence} valid cells where "
                "confidence is zero"
            )
        band_reports.append(
            {
                "name": name,
                "valid_cells": valid_count,
                "nodata_cells": missing_count,
                "missing_with_positive_confidence": (
                    missing_with_positive_confidence
                    if name != "confidence"
                    else None
                ),
                "valid_with_zero_confidence": (
                    valid_with_zero_confidence if name != "confidence" else None
                ),
            }
        )

    output.FlushCache()
    output = None
    remove_own_output(OUTPUT_PATH)
    TEMP_PATH.replace(OUTPUT_PATH)
    assert_matches_grid(OUTPUT_PATH)
    return band_reports


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    source_before = {
        "size": SOURCE_PATH.stat().st_size,
        "mtime_ns": SOURCE_PATH.stat().st_mtime_ns,
    }
    source = gdal.OpenEx(str(SOURCE_PATH), gdal.OF_RASTER | gdal.OF_READONLY)
    if source is None:
        raise RuntimeError(f"cannot open {SOURCE_PATH} read-only")

    discarded, confidence = verify_source(source)
    band_reports = create_output(source, confidence)
    source = None
    source_after = {
        "size": SOURCE_PATH.stat().st_size,
        "mtime_ns": SOURCE_PATH.stat().st_mtime_ns,
    }
    if source_before != source_after:
        raise RuntimeError("source size or modification time changed")

    report = {
        "status": "passed",
        "operation": "window-copy columns 0:1812; no interpolation",
        "source": str(SOURCE_PATH.relative_to(ROOT)),
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
        "source_unchanged": True,
        "source_shape_rows_cols": [ROWS, COLS + 1],
        "output_shape_rows_cols": [ROWS, COLS],
        "physical_dtype": "float32",
        "physical_nodata": NODATA,
        "discarded_eastern_column": discarded,
        "bands": band_reports,
        "missingness_policy": (
            "nutrient-specific nulls preserved independently of confidence"
        ),
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
