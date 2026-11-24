#!/usr/bin/env python3
"""Stage the owner-approved AfricanSoilGrids phosphorus without alteration.

The owner confirmed mg/100kg on 2026-11-24. Preserve those values exactly:
no additional scaling, back-transformation, resampling, clipping, or imputation.
The owner subsequently identified AfricanSoilGrids (not global SoilGrids) and
explicitly approved using the supplied values as-is, with no scaling or other
transform. Only source gaps and the common Ethiopia mask affect stored pixels.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config.grid import (  # noqa: E402
    COLS,
    CRS,
    NODATA,
    ROWS,
    TRANSFORM,
    assert_matches_grid,
)

SOURCE_PATH = ROOT / "data/inputs/SoilGrids/p_0-30cm_mean_1km_4326.tif"
MASK_PATH = ROOT / "data/interim/admin_grid/ethiopia_mask.tif"
HERE = ROOT / "data/interim/phosphorus"
OUTPUT_PATH = HERE / "phosphorus_staged.tif"
REPORT_PATH = HERE / "phosphorus_report.json"
OUTPUT_TEMP_PATH = HERE / "phosphorus_staged.tmp.tif"
REPORT_TEMP_PATH = HERE / "phosphorus_report.tmp.json"
NAME = "sg_p"
UNIT = "mg/100kg"
SOURCE_LABEL = "AfricanSoilGrids (owner-confirmed; not global SoilGrids)"
OWNER_APPROVAL = {
    "date": "2026-11-24",
    "authority": "owner confirmation in implementation conversation",
    "product": "AfricanSoilGrids",
    "unit": UNIT,
    "value_policy": "Use supplied values as-is; no scaling or change required.",
    "required_transform": "identity; no scaling or back-transformation",
}

gdal.UseExceptions()


def relative(path: Path) -> str:
    """Keep project paths portable, while permitting temporary test fixtures."""
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def open_raster(path: Path):
    dataset = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
    if dataset is None:
        raise RuntimeError(f"cannot open raster read-only: {relative(path)}")
    return dataset


def band_by_name(dataset, name: str, *, unnamed_alias: str | None = None):
    """Resolve descriptions; alias an unnamed singleton only when explicit.

    The inspected phosphorus source has no description. Its verified filename
    binds its sole band to sg_p; there is no assumed multiband source position.
    Named masks and outputs always resolve their actual descriptions.
    """
    if unnamed_alias is not None and dataset.RasterCount != 1:
        raise ValueError("an unnamed source alias requires exactly one band")
    bands = {}
    for index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(index)
        description = band.GetDescription() or unnamed_alias
        if not description or description in bands:
            raise ValueError(f"missing or duplicate raster band name: {description!r}")
        bands[description] = band
    if name not in bands:
        raise ValueError(f"missing raster band {name!r}; found {list(bands)!r}")
    return bands[name]


def fingerprints() -> dict:
    """Hash source/mask files and sidecars without updating their metadata."""
    paths = set()
    for path in (SOURCE_PATH, MASK_PATH):
        dataset = open_raster(path)
        paths.update(Path(filename) for filename in dataset.GetFileList())
        dataset = None
    records = {}
    for path in sorted(paths):
        with path.open("rb") as source:
            checksum = hashlib.file_digest(source, "sha256").hexdigest()
        records[relative(path)] = {
            "size_bytes": path.stat().st_size,
            "sha256": checksum,
        }
    return records


def read_source() -> tuple[np.ndarray, dict]:
    assert_matches_grid(SOURCE_PATH)
    dataset = open_raster(SOURCE_PATH)
    if dataset.RasterCount != 1:
        raise ValueError("phosphorus source must contain exactly one band")
    band = band_by_name(dataset, NAME, unnamed_alias=NAME)
    if band.DataType != gdal.GDT_Float32:
        raise ValueError("phosphorus source must be Float32")
    if band.GetNoDataValue() != NODATA:
        raise ValueError(f"phosphorus source nodata must equal canonical {NODATA}")
    if band.GetScale() not in (None, 1.0) or band.GetOffset() not in (None, 0.0):
        raise ValueError("phosphorus has an embedded scale/offset; review before copying")
    for unit in (band.GetUnitType(), band.GetMetadataItem("unit")):
        if unit and unit != UNIT:
            raise ValueError(f"embedded phosphorus unit {unit!r} disagrees with {UNIT!r}")
    values = band.ReadAsArray()
    if values is None or values.shape != (ROWS, COLS):
        raise RuntimeError("failed reading phosphorus on the canonical grid")
    metadata = {
        "file": relative(SOURCE_PATH),
        "band_description": band.GetDescription(),
        "band_binding": "verified singleton filename alias to sg_p",
        "embedded_unit": band.GetUnitType() or None,
        "embedded_scale": band.GetScale(),
        "embedded_offset": band.GetOffset(),
        "nodata": band.GetNoDataValue(),
        "upstream_product": "AfricanSoilGrids",
        "required_scaling_and_back_transform": "none; owner explicitly approved supplied values as-is",
        "provenance_authority": OWNER_APPROVAL,
        "product_version": None,
        "product_version_note": "Not supplied; no version or processing history invented.",
    }
    band = None
    dataset = None
    return values, metadata


def read_mask() -> np.ndarray:
    assert_matches_grid(MASK_PATH)
    dataset = open_raster(MASK_PATH)
    if dataset.RasterCount != 1:
        raise ValueError("ethiopia_mask must contain exactly one band")
    values = band_by_name(dataset, "ethiopia_mask").ReadAsArray()
    dataset = None
    if values is None or values.shape != (ROWS, COLS):
        raise RuntimeError("failed reading the country mask")
    if not np.all(np.isin(values, (0, 1))):
        raise ValueError("ethiopia_mask must contain only 0 and 1")
    return values == 1


def masked_copy(values: np.ndarray, inside: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Preserve every valid value, including zero, and every source gap."""
    if values.ndim != 2 or inside.shape != values.shape or inside.dtype != np.bool_:
        raise ValueError("source and boolean country mask must have matching 2-D shapes")
    valid = np.isfinite(values) & (values != NODATA)
    staged = np.full(values.shape, NODATA, dtype=np.float32)
    staged[inside & valid] = values[inside & valid]
    return staged, valid


def write_output(values: np.ndarray) -> None:
    output = gdal.GetDriverByName("GTiff").Create(
        str(OUTPUT_TEMP_PATH), COLS, ROWS, 1, gdal.GDT_Float32,
        options=[
            "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE", "PREDICTOR=3", "INTERLEAVE=BAND",
        ],
    )
    if output is None:
        raise RuntimeError(f"cannot create {relative(OUTPUT_TEMP_PATH)}")
    output.SetGeoTransform(TRANSFORM.to_gdal())
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(CRS)
    output.SetSpatialRef(spatial_ref)
    output.SetMetadataItem("AREA_OR_POINT", "Area")
    output.SetMetadataItem("role", "owner-approved AfricanSoilGrids phosphorus staging")
    output.SetMetadataItem("stack_eligible", "true")
    output.SetMetadataItem("publication_ready", "true")
    # This alias is only used while creating the new unnamed singleton. Reads
    # below require the persisted description rather than assuming its index.
    band = band_by_name(output, NAME, unnamed_alias=NAME)
    band.SetDescription(NAME)
    band.SetNoDataValue(NODATA)
    band.SetUnitType(UNIT)
    for key, value in {
        "unit": UNIT,
        "source": SOURCE_LABEL,
        "source_file": relative(SOURCE_PATH),
        "logical_dtype": "float32",
        "resolution": "1 km common-grid cell; upstream native resolution unverified",
        "resample": "none",
        "transform_applied": "none; retain owner-confirmed mg/100kg",
        "provenance_status": "owner-approved AfricanSoilGrids; supplied values require no transform",
    }.items():
        band.SetMetadataItem(key, value)
    band.WriteArray(values)
    output.FlushCache()
    band = None
    output = None


def validate_output(expected: np.ndarray, inside: np.ndarray) -> dict:
    assert_matches_grid(OUTPUT_TEMP_PATH)
    dataset = open_raster(OUTPUT_TEMP_PATH)
    if dataset.RasterCount != 1:
        raise ValueError("staged phosphorus must contain exactly one band")
    band = band_by_name(dataset, NAME)
    if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != NODATA:
        raise ValueError("staged phosphorus storage contract changed")
    if band.GetUnitType() != UNIT or band.GetMetadataItem("unit") != UNIT:
        raise ValueError("staged phosphorus unit did not round-trip")
    if dataset.GetMetadataItem("stack_eligible") != "true":
        raise ValueError("approved phosphorus lost its recorded inclusion approval")
    observed = band.ReadAsArray()
    if observed is None or observed.shape != expected.shape:
        raise RuntimeError("failed reopening staged phosphorus")
    mismatches = int(np.count_nonzero(observed.view(np.uint32) != expected.view(np.uint32)))
    outside_populated = int(np.count_nonzero((observed != NODATA) & ~inside))
    if mismatches or outside_populated:
        raise ValueError(
            f"phosphorus copy mismatches={mismatches}; populated outside={outside_populated}"
        )
    band = None
    dataset = None
    return {"bitwise_round_trip_mismatches": mismatches,
            "populated_outside_cells_after_mask": outside_populated}


def main() -> None:
    """Overwrite only this step's artifacts after all technical checks pass."""
    OUTPUT_TEMP_PATH.unlink(missing_ok=True)
    REPORT_TEMP_PATH.unlink(missing_ok=True)
    try:
        before = fingerprints()
        values, source_metadata = read_source()
        inside = read_mask()
        staged, valid = masked_copy(values, inside)
        selected = values[inside & valid].astype(np.float64)
        if not selected.size:
            raise ValueError("phosphorus has no valid cells inside Ethiopia")
        write_output(staged)
        checks = validate_output(staged, inside)
        if fingerprints() != before:
            raise RuntimeError("source or mask files changed during phosphorus staging")

        statistics = {
            "valid_inside_cells": int(selected.size),
            "missing_inside_cells": int(np.count_nonzero(inside & ~valid)),
            "nonfinite_inside_cells": int(np.count_nonzero(inside & ~np.isfinite(values))),
            "populated_outside_cells_before_mask": int(np.count_nonzero(~inside & valid)),
            "min": float(selected.min()),
            "median": float(np.median(selected)),
            "mean": float(selected.mean()),
            "max": float(selected.max()),
        }
        report = {
            "status": "passed",
            "technical_validation_status": "passed",
            "stack_eligible": True,
            "publication_ready": True,
            "output": relative(OUTPUT_PATH),
            "grid_check": "passed",
            "source": source_metadata,
            "source_fingerprints": before,
            "sources_unchanged": True,
            "unit_confirmation": {
                "unit": UNIT,
                "authority": "owner confirmation in implementation conversation",
                "date": "2026-11-24",
            },
            "owner_approval": OWNER_APPROVAL,
            "copy_policy": "no resampling, scaling, back-transformation, clipping, or imputation",
            "mask": {
                "path": relative(MASK_PATH),
                "inside_cells": int(np.count_nonzero(inside)),
                "outside_cells": int(np.count_nonzero(~inside)),
            },
            "bands": [{
                "name": NAME,
                "group": "soil",
                "description": "AfricanSoilGrids phosphorus; supplied 0-30 cm common-grid product",
                "unit": UNIT,
                "physical_dtype": "float32",
                "logical_dtype": "float32",
                "nodata": NODATA,
                "source": SOURCE_LABEL,
                "source_file": relative(SOURCE_PATH),
                "resample": "none; exact-grid masked copy",
                "transform_applied": "none",
                "overview_resample": "average",
                "statistics": statistics,
                **checks,
            }],
            # These are the existing specification's diagnostic limits, not
            # owner-approved scientific bounds or scoring thresholds.
            "provisional_spec_diagnostic": {
                "source": "Specification/01_grid_and_bands.md section 4",
                "unit": UNIT,
                "range": [70, 4000],
                "below_range_cells": int(np.count_nonzero(selected < 70)),
                "above_range_cells": int(np.count_nonzero(selected > 4000)),
                "policy": "report only; no clipping or scientific acceptance implied",
            },
            "dimensional_reference_only": {
                "note": "Dividing mg/100kg by 100 gives mg/kg; the raster is NOT converted.",
                "unit": "mg/kg",
                "min": statistics["min"] / 100,
                "median": statistics["median"] / 100,
                "max": statistics["max"] / 100,
            },
            "publication_blockers": [],
        }
        text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        REPORT_TEMP_PATH.write_text(text)
        OUTPUT_TEMP_PATH.replace(OUTPUT_PATH)
        REPORT_TEMP_PATH.replace(REPORT_PATH)
        print(text, end="")
    finally:
        OUTPUT_TEMP_PATH.unlink(missing_ok=True)
        REPORT_TEMP_PATH.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
