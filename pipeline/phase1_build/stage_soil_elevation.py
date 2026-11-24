#!/usr/bin/env python3
"""Stage approved SoilGrids and elevation bands without resampling.

Every source must already match the canonical grid. Values are copied exactly
inside Ethiopia, source missing values remain missing, and all values outside
the approved country mask become the physical -9999 sentinel. Phosphorus is
deliberately excluded until its product-specific provenance is approved.
"""

from __future__ import annotations

import json
import math
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

SOIL_DIR = ROOT / "data/inputs/SoilGrids"
ELEVATION_DIR = ROOT / "data/inputs/Elevation"
MASK_PATH = ROOT / "data/interim/admin_grid/ethiopia_mask.tif"
PHOSPHORUS_PATH = SOIL_DIR / "p_0-30cm_mean_1km_4326.tif"
HERE = ROOT / "data/interim/soil_elevation"
OUTPUT_PATH = HERE / "soil_elevation_staged.tif"
REPORT_PATH = HERE / "soil_elevation_report.json"
OUTPUT_TEMP_PATH = HERE / "soil_elevation_staged.tmp.tif"
REPORT_TEMP_PATH = HERE / "soil_elevation_report.tmp.json"

BANDS = (
    {
        "name": "sg_phh2o",
        "description": "Soil pH in water, 0-30 cm depth-weighted mean",
        "unit": "pH",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "phh2o_0-30cm_mean_1km_4326.tif",
        "valid_min": 3.0,
        "valid_max": 10.0,
        "notes": "Already converted to application units upstream.",
    },
    {
        "name": "sg_sand",
        "description": "Soil sand content, 0-30 cm depth-weighted mean",
        "unit": "%",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "sand_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 100.0,
        "notes": "Already converted to percent upstream.",
    },
    {
        "name": "sg_clay",
        "description": "Soil clay content, 0-30 cm depth-weighted mean",
        "unit": "%",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "clay_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 100.0,
        "notes": "Already converted to percent upstream.",
    },
    {
        "name": "sg_silt",
        "description": "Soil silt content, 0-30 cm depth-weighted mean",
        "unit": "%",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "silt_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 100.0,
        "notes": "Already converted to percent upstream.",
    },
    {
        "name": "sg_cec",
        "description": "Soil cation exchange capacity, 0-30 cm depth-weighted mean",
        "unit": "cmol(c)/kg",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "cec_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 100.0,
        "notes": "Already converted to application units upstream.",
    },
    {
        "name": "sg_soc",
        "description": "Soil organic carbon, 0-30 cm depth-weighted mean",
        "unit": "g/kg",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "soc_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 250.0,
        "notes": "Available application organic-carbon input; no es_oc exists.",
    },
    {
        "name": "sg_n",
        "description": "Soil total nitrogen, 0-30 cm depth-weighted mean",
        "unit": "mg/kg",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "n_0-30cm_mean_1km_4326.tif",
        "valid_min": 500.0,
        "valid_max": 5000.0,
        "notes": "Sole application nitrogen input; already mg/kg; no scaling applied.",
    },
    {
        "name": "sg_cfvo",
        "description": "Soil coarse fragments, 0-30 cm depth-weighted mean",
        "unit": "vol %",
        "group": "soil",
        "source": "SoilGrids 2.0 (ISRIC)",
        "path": SOIL_DIR / "cfvo_0-30cm_mean_1km_4326.tif",
        "valid_min": 0.0,
        "valid_max": 100.0,
        "notes": "Already converted to volume percent upstream.",
    },
    {
        "name": "elevation_m",
        "description": "Mean surface elevation on the 1 km common grid",
        "unit": "m",
        "group": "terrain",
        "source": "Copernicus DEM GLO-90",
        "path": ELEVATION_DIR / "ethiopia_Copernicus-DEM_1KM.tif",
        "valid_min": -200.0,
        "valid_max": 4700.0,
        "notes": "EGM2008 geoid; surface model including canopy; 1 km aggregate.",
    },
)

gdal.UseExceptions()


def relative(path: Path) -> str:
    """Return a repository-relative path for metadata and diagnostics."""
    return path.relative_to(ROOT).as_posix()


def open_raster(path: Path):
    """Open one source raster read-only."""
    dataset = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
    if dataset is None:
        raise RuntimeError(f"cannot open raster read-only: {relative(path)}")
    return dataset


def validate_input_inventory() -> None:
    """Fail rather than silently omit an unexpected SoilGrids raster."""
    selected_soil = {
        band["path"] for band in BANDS if band["group"] == "soil"
    }
    expected_soil = selected_soil | {PHOSPHORUS_PATH}
    observed_soil = set(SOIL_DIR.glob("*.tif"))
    if observed_soil != expected_soil:
        missing = sorted(relative(path) for path in expected_soil - observed_soil)
        unexpected = sorted(
            relative(path) for path in observed_soil - expected_soil
        )
        raise ValueError(
            "SoilGrids inventory changed; review before staging: "
            f"missing={missing}, unexpected={unexpected}"
        )
    if not PHOSPHORUS_PATH.is_file():
        raise FileNotFoundError(relative(PHOSPHORUS_PATH))
    if any(band["name"] == "sg_p" for band in BANDS):
        raise AssertionError("gated phosphorus must not enter this staging build")


def read_mask() -> np.ndarray:
    """Read and validate the previously approved country mask."""
    assert_matches_grid(MASK_PATH)
    dataset = open_raster(MASK_PATH)
    if dataset.RasterCount != 1:
        raise ValueError(f"{relative(MASK_PATH)} must contain one band")
    mask = dataset.GetRasterBand(1).ReadAsArray()
    dataset = None
    if mask is None or mask.shape != (ROWS, COLS):
        raise RuntimeError(f"failed reading {relative(MASK_PATH)}")
    if not np.all(np.isin(mask, (0, 1))):
        raise ValueError("ethiopia_mask contains values other than 0 and 1")
    return mask == 1


def create_output():
    """Create the temporary uniform-float32 staging dataset."""
    output = gdal.GetDriverByName("GTiff").Create(
        str(OUTPUT_TEMP_PATH),
        COLS,
        ROWS,
        len(BANDS),
        gdal.GDT_Float32,
        options=[
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
            "INTERLEAVE=BAND",
            "NUM_THREADS=ALL_CPUS",
        ],
    )
    if output is None:
        raise RuntimeError(f"cannot create {relative(OUTPUT_TEMP_PATH)}")
    output.SetGeoTransform(TRANSFORM.to_gdal())
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(CRS)
    output.SetSpatialRef(spatial_ref)
    output.SetMetadataItem("AREA_OR_POINT", "Area")
    output.SetMetadataItem("role", "Phase 1 soil and elevation staging")
    output.SetMetadataItem("resample", "none; exact-grid copies")
    output.SetMetadataItem("phosphorus", "excluded pending provenance approval")
    return output


def source_valid_mask(values: np.ndarray, nodata: float) -> np.ndarray:
    """Return finite values not equal to a source's declared nodata."""
    if math.isnan(nodata):
        tagged_missing = np.isnan(values)
    else:
        tagged_missing = values == nodata
    return np.isfinite(values) & ~tagged_missing


def stage_band(output, index: int, definition: dict, inside: np.ndarray) -> dict:
    """Copy one exact-grid source band, mask it, and return an audit record."""
    source_path = definition["path"]
    assert_matches_grid(source_path)
    source = open_raster(source_path)
    if source.RasterCount != 1:
        raise ValueError(
            f"{relative(source_path)} has {source.RasterCount} bands; expected 1"
        )
    source_band = source.GetRasterBand(1)
    if gdal.GetDataTypeName(source_band.DataType) != "Float32":
        raise ValueError(f"{relative(source_path)} must be Float32")
    source_nodata = source_band.GetNoDataValue()
    if source_nodata is None or not math.isclose(
        source_nodata, NODATA, rel_tol=0, abs_tol=0
    ):
        raise ValueError(
            f"{relative(source_path)} nodata {source_nodata!r} != {NODATA}"
        )
    source_scale = source_band.GetScale()
    source_offset = source_band.GetOffset()
    if source_scale not in (None, 1.0) or source_offset not in (None, 0.0):
        raise ValueError(
            f"{relative(source_path)} has an unapplied scale/offset: "
            f"{source_scale!r}/{source_offset!r}"
        )

    values = source_band.ReadAsArray()
    if values is None or values.shape != (ROWS, COLS):
        raise RuntimeError(f"failed reading {relative(source_path)}")
    valid = source_valid_mask(values, source_nodata)
    valid_inside = valid & inside
    staged = np.full((ROWS, COLS), NODATA, dtype=np.float32)
    staged[valid_inside] = values[valid_inside]

    valid_values = values[valid_inside].astype(np.float64, copy=False)
    if valid_values.size == 0:
        raise ValueError(f"{definition['name']} has no valid cells inside Ethiopia")

    output_band = output.GetRasterBand(index)
    output_band.SetNoDataValue(NODATA)
    output_band.SetDescription(definition["name"])
    output_band.SetMetadataItem("unit", definition["unit"])
    output_band.SetMetadataItem("source", definition["source"])
    output_band.SetMetadataItem("source_file", relative(source_path))
    output_band.SetMetadataItem("resample", "none")
    output_band.SetMetadataItem("logical_dtype", "float32")
    output_band.WriteArray(staged)

    report = {
        "index": index,
        "name": definition["name"],
        "description": definition["description"],
        "group": definition["group"],
        "unit": definition["unit"],
        "physical_dtype": "float32",
        "logical_dtype": "float32",
        "nodata": NODATA,
        "source": definition["source"],
        "source_file": relative(source_path),
        "source_band_description": source_band.GetDescription(),
        "source_unit_tag": source_band.GetUnitType(),
        "source_scale": source_scale,
        "source_offset": source_offset,
        "resample": "none; exact-grid copy",
        "valid_min": definition["valid_min"],
        "valid_max": definition["valid_max"],
        "notes": definition["notes"],
        "statistics": {
            "valid_inside_cells": int(valid_values.size),
            "missing_inside_cells": int(np.count_nonzero(inside & ~valid)),
            "nonfinite_inside_cells": int(
                np.count_nonzero(inside & ~np.isfinite(values))
            ),
            "populated_outside_cells_before_mask": int(
                np.count_nonzero(~inside & valid)
            ),
            "min": float(valid_values.min()),
            "median": float(np.median(valid_values)),
            "mean": float(valid_values.mean(dtype=np.float64)),
            "max": float(valid_values.max()),
            "below_sanity_range_cells": int(
                np.count_nonzero(valid_values < definition["valid_min"])
            ),
            "above_sanity_range_cells": int(
                np.count_nonzero(valid_values > definition["valid_max"])
            ),
        },
    }
    source = None
    return report


def validate_round_trip(inside: np.ndarray, reports: list[dict]) -> None:
    """Prove staged values equal each source exactly wherever both are valid."""
    assert_matches_grid(OUTPUT_TEMP_PATH)
    staged = open_raster(OUTPUT_TEMP_PATH)
    if staged.RasterCount != len(BANDS):
        raise ValueError(
            f"staged band count {staged.RasterCount} != {len(BANDS)}"
        )

    descriptions = []
    for index, (definition, report) in enumerate(zip(BANDS, reports), start=1):
        source = open_raster(definition["path"])
        source_band = source.GetRasterBand(1)
        source_values = source_band.ReadAsArray()
        valid = source_valid_mask(source_values, source_band.GetNoDataValue())
        expected = np.full((ROWS, COLS), NODATA, dtype=np.float32)
        expected[valid & inside] = source_values[valid & inside]

        staged_band = staged.GetRasterBand(index)
        staged_values = staged_band.ReadAsArray()
        mismatches = int(np.count_nonzero(staged_values != expected))
        report["round_trip_mismatches"] = mismatches
        if mismatches:
            raise ValueError(
                f"{definition['name']} has {mismatches} source-copy mismatches"
            )
        descriptions.append(staged_band.GetDescription())
        source = None

    expected_descriptions = [band["name"] for band in BANDS]
    if descriptions != expected_descriptions:
        raise ValueError(
            f"staged descriptions {descriptions} != {expected_descriptions}"
        )
    if "sg_p" in descriptions:
        raise AssertionError("phosphorus entered the staged raster")
    staged = None


def main() -> None:
    """Build and atomically publish the staged raster and report."""
    OUTPUT_TEMP_PATH.unlink(missing_ok=True)
    REPORT_TEMP_PATH.unlink(missing_ok=True)
    try:
        validate_input_inventory()
        inside = read_mask()
        output = create_output()
        reports = [
            stage_band(output, index, definition, inside)
            for index, definition in enumerate(BANDS, start=1)
        ]
        output.FlushCache()
        output = None

        validate_round_trip(inside, reports)
        report = {
            "status": "passed",
            "output": relative(OUTPUT_PATH),
            "grid_check": "passed",
            "physical_dtype": "float32",
            "nodata": NODATA,
            "mask": {
                "path": relative(MASK_PATH),
                "inside_cells": int(np.count_nonzero(inside)),
                "outside_cells": int(np.count_nonzero(~inside)),
                "outside_values_after_staging": "all nodata",
            },
            "copy_policy": "no resampling, scaling, clipping, or imputation",
            "band_order": [band["name"] for band in BANDS],
            "bands": reports,
            "excluded": [
                {
                    "name": "sg_p",
                    "source_file": relative(PHOSPHORUS_PATH),
                    "reason": (
                        "blocked pending exact product, unit, scaling, and "
                        "back-transformation approval"
                    ),
                }
            ],
            "sources_unchanged": True,
        }
        REPORT_TEMP_PATH.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        OUTPUT_TEMP_PATH.replace(OUTPUT_PATH)
        assert_matches_grid(OUTPUT_PATH)
        REPORT_TEMP_PATH.replace(REPORT_PATH)
        print(json.dumps(report, indent=2, sort_keys=True))
    finally:
        OUTPUT_TEMP_PATH.unlink(missing_ok=True)
        REPORT_TEMP_PATH.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
