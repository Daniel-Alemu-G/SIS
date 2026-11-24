#!/usr/bin/env python3
"""Stage land-cover, seasonal rainfall, and temperature bands.

The supplied 1 km products already match the canonical grid, so this step
never resamples them. It preserves legal land-cover class codes, validates the
rainfall companion metadata, maps temperature descriptions explicitly to the
canonical ``_c`` names, and applies the approved Ethiopia mask.
"""

from __future__ import annotations

import hashlib
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

LANDCOVER_PATH = ROOT / "data/inputs/Landcover/WorldCover_Ethiopia_1km.tif"
RAINFALL_PATH = (
    ROOT / "data/inputs/Rainfall/seasonal/rainfall_seasonal_1km.tif"
)
RAINFALL_METADATA_PATH = (
    ROOT / "data/inputs/Rainfall/seasonal/rainfall_seasonal_1km.json"
)
TEMPERATURE_PATH = ROOT / "data/inputs/temprature/tmean_seasonal.tif"
MASK_PATH = ROOT / "data/interim/admin_grid/ethiopia_mask.tif"
HERE = ROOT / "data/interim/landcover_climate"
OUTPUT_PATH = HERE / "landcover_climate_staged.tif"
REPORT_PATH = HERE / "landcover_climate_report.json"
OUTPUT_TEMP_PATH = HERE / "landcover_climate_staged.tmp.tif"
REPORT_TEMP_PATH = HERE / "landcover_climate_report.tmp.json"

LANDCOVER_CLASSES = {
    10: "tree cover",
    20: "shrubland",
    30: "grassland",
    40: "cropland",
    50: "built-up",
    60: "bare or sparse vegetation",
    70: "snow and ice",
    80: "permanent water bodies",
    90: "herbaceous wetland",
    95: "mangroves",
    100: "moss and lichen",
}

RAINFALL_NAMES = (
    "kiremt_mean_mm",
    "kiremt_std_mm",
    "kiremt_cv",
    "belg_mean_mm",
    "belg_std_mm",
    "belg_cv",
    "bega_mean_mm",
    "bega_std_mm",
    "bega_cv",
)

SEASON_MONTHS = {
    "kiremt": [6, 7, 8, 9],
    "belg": [2, 3, 4, 5],
    "bega": [10, 11, 12, 1],
}

RAINFALL_RANGES = {
    "kiremt_mean_mm": (0.0, 2500.0),
    "kiremt_std_mm": (0.0, None),
    "kiremt_cv": (0.0, 2.0),
    "belg_mean_mm": (0.0, 1000.0),
    "belg_std_mm": (0.0, None),
    "belg_cv": (0.0, 2.0),
    "bega_mean_mm": (0.0, 800.0),
    "bega_std_mm": (0.0, None),
    "bega_cv": (0.0, 2.0),
}

TEMPERATURE_MAPPINGS = (
    ("tmean_kiremt_c", "tmean_kiremt", "Kiremt mean temperature"),
    ("tmean_bega_c", "tmean_bega", "Bega mean temperature"),
    ("tmean_belg_c", "tmean_belg", "Belg mean temperature"),
)

gdal.UseExceptions()


def relative(path: Path) -> str:
    """Return a stable repository-relative path."""
    return path.relative_to(ROOT).as_posix()


def open_raster(path: Path):
    """Open a raster read-only or fail with a clear path."""
    dataset = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
    if dataset is None:
        raise RuntimeError(f"cannot open raster read-only: {relative(path)}")
    return dataset


def band_indices_by_name(dataset) -> dict[str, int]:
    """Resolve band indices from descriptions, rejecting ambiguous names."""
    indices = {}
    for index in range(1, dataset.RasterCount + 1):
        name = dataset.GetRasterBand(index).GetDescription()
        if not name or name in indices:
            raise ValueError(f"missing or duplicate source band name: {name!r}")
        indices[name] = index
    return indices


def source_fingerprints() -> dict:
    """Fingerprint files and sidecars without writing to input directories."""
    paths = {RAINFALL_METADATA_PATH, MASK_PATH}
    for path in (LANDCOVER_PATH, RAINFALL_PATH, TEMPERATURE_PATH):
        dataset = open_raster(path)
        paths.update(Path(filename) for filename in dataset.GetFileList())
        dataset = None
    fingerprints = {}
    for path in sorted(paths):
        checksum = hashlib.sha256()
        with path.open("rb") as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                checksum.update(chunk)
        fingerprints[relative(path)] = {
            "size_bytes": path.stat().st_size,
            "sha256": checksum.hexdigest(),
        }
    return fingerprints


def validate_source_layouts() -> None:
    """Validate complete source band layouts before reading any values."""
    for path in (LANDCOVER_PATH, RAINFALL_PATH, TEMPERATURE_PATH):
        assert_matches_grid(path)

    landcover = open_raster(LANDCOVER_PATH)
    if landcover.RasterCount != 1:
        raise ValueError(f"{relative(LANDCOVER_PATH)} must contain one band")
    if landcover.GetRasterBand(1).GetDescription() != "Map":
        raise ValueError("land-cover source band description must be 'Map'")
    landcover = None

    rainfall = open_raster(RAINFALL_PATH)
    rainfall_descriptions = [
        rainfall.GetRasterBand(index).GetDescription()
        for index in range(1, rainfall.RasterCount + 1)
    ]
    if rainfall_descriptions != list(RAINFALL_NAMES):
        raise ValueError(
            f"rainfall descriptions {rainfall_descriptions} != "
            f"{list(RAINFALL_NAMES)}"
        )
    rainfall = None

    temperature = open_raster(TEMPERATURE_PATH)
    temperature_descriptions = [
        temperature.GetRasterBand(index).GetDescription()
        for index in range(1, temperature.RasterCount + 1)
    ]
    expected_temperature = ["tmean_kiremt", "tmean_belg", "tmean_bega"]
    if temperature_descriptions != expected_temperature:
        raise ValueError(
            f"temperature descriptions {temperature_descriptions} != "
            f"{expected_temperature}"
        )
    temperature = None


def load_rainfall_metadata() -> list[dict]:
    """Validate and return the supplied seasonal-rainfall provenance records."""
    records = json.loads(RAINFALL_METADATA_PATH.read_text())
    if not isinstance(records, list) or len(records) != len(RAINFALL_NAMES):
        raise ValueError("rainfall companion JSON must contain nine records")

    for expected_index, (expected_name, record) in enumerate(
        zip(RAINFALL_NAMES, records), start=1
    ):
        if not isinstance(record, dict):
            raise ValueError(f"rainfall metadata record {expected_index} is not an object")
        if record.get("band") != expected_index:
            raise ValueError(
                f"rainfall metadata band {record.get('band')!r} != "
                f"{expected_index}"
            )
        if record.get("name") != expected_name:
            raise ValueError(
                f"rainfall metadata name {record.get('name')!r} != "
                f"{expected_name!r}"
            )
        season = record.get("season")
        if (
            season != expected_name.split("_", 1)[0]
            or record.get("months") != SEASON_MONTHS.get(season)
        ):
            raise ValueError(
                f"{expected_name} months {record.get('months')!r} do not "
                f"match the approved {season!r} season"
            )
        years = record.get("years_used")
        if (
            not isinstance(years, list)
            or not all(type(year) is int for year in years)
            or len(years) != len(set(years))
            or years != sorted(years)
            or record.get("n_years") != len(years)
            or record.get("n_years") != 15
        ):
            raise ValueError(f"{expected_name} has inconsistent year metadata")

        expected_stat = (
            "mean_mm"
            if expected_name.endswith("_mean_mm")
            else "std_mm"
            if expected_name.endswith("_std_mm")
            else "cv"
        )
        if record.get("stat") != expected_stat:
            raise ValueError(
                f"{expected_name} stat {record.get('stat')!r} != "
                f"{expected_stat!r}"
            )
    for season in SEASON_MONTHS:
        years_by_stat = [
            record["years_used"] for record in records if record["season"] == season
        ]
        if any(years != years_by_stat[0] for years in years_by_stat[1:]):
            raise ValueError(f"{season} rainfall statistics use different season-years")
    return records


def build_band_definitions(rainfall_records: list[dict]) -> list[dict]:
    """Create manifest-ready definitions in the fixed final group order."""
    source_indices = {}
    for path in (LANDCOVER_PATH, RAINFALL_PATH, TEMPERATURE_PATH):
        dataset = open_raster(path)
        source_indices[path] = band_indices_by_name(dataset)
        dataset = None
    definitions = [
        {
            "name": "landcover",
            "description": "Dominant ESA WorldCover class on the 1 km grid",
            "group": "landcover",
            "unit": "class code",
            "logical_dtype": "uint8",
            "source": "ESA WorldCover",
            "path": LANDCOVER_PATH,
            "source_band": source_indices[LANDCOVER_PATH]["Map"],
            "source_description": "Map",
            "source_dtype": "Byte",
            "source_nodata": 0.0,
            "valid_min": 10.0,
            "valid_max": 100.0,
            "overview_resample": "mode",
            "notes": (
                "Dominant class only; preserve integer codes; never average or "
                "treat as cropland percentage."
            ),
            "provenance": None,
        }
    ]

    for record in rainfall_records:
        name = record["name"]
        valid_min, valid_max = RAINFALL_RANGES[name]
        unit = "fraction" if name.endswith("_cv") else "mm"
        stat_description = (
            "coefficient of variation"
            if name.endswith("_cv")
            else "sample standard deviation"
            if name.endswith("_std_mm")
            else "mean"
        )
        definitions.append(
            {
                "name": name,
                "description": (
                    f"{record['season'].capitalize()} seasonal rainfall "
                    f"{stat_description}"
                ),
                "group": "rain",
                "unit": unit,
                "logical_dtype": "float32",
                "source": "CHIRPS seasonal aggregate (supplied product)",
                "path": RAINFALL_PATH,
                "source_band": source_indices[RAINFALL_PATH][name],
                "source_description": name,
                "source_dtype": "Float32",
                "source_nodata": NODATA,
                "valid_min": valid_min,
                "valid_max": valid_max,
                "overview_resample": "average",
                "notes": (
                    "Precomputed seasonal product; no monthly rainfall rebuild."
                ),
                "provenance": record,
            }
        )

    for canonical_name, source_name, description in (
        TEMPERATURE_MAPPINGS
    ):
        definitions.append(
            {
                "name": canonical_name,
                "description": description,
                "group": "temperature",
                "unit": "degrees Celsius",
                "logical_dtype": "float32",
                "source": "Seasonal mean temperature (supplied product)",
                "path": TEMPERATURE_PATH,
                "source_band": source_indices[TEMPERATURE_PATH][source_name],
                "source_description": source_name,
                "source_dtype": "Float32",
                "source_nodata": NODATA,
                "valid_min": -20.0,
                "valid_max": 70.0,
                "overview_resample": "average",
                "notes": "Mapped explicitly from the source description.",
                "provenance": None,
            }
        )
    return definitions


def read_mask() -> np.ndarray:
    """Read and validate the approved Ethiopia mask."""
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


def create_output(band_count: int):
    """Create the temporary uniform-float32 staging dataset."""
    output = gdal.GetDriverByName("GTiff").Create(
        str(OUTPUT_TEMP_PATH),
        COLS,
        ROWS,
        band_count,
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
    output.SetMetadataItem("role", "Phase 1 land-cover and climate staging")
    output.SetMetadataItem("resample", "none; exact-grid copies")
    return output


def valid_mask(values: np.ndarray, nodata: float) -> np.ndarray:
    """Return finite values that are not the source nodata value."""
    tagged_missing = np.isnan(values) if math.isnan(nodata) else values == nodata
    return np.isfinite(values) & ~tagged_missing


def range_count(values: np.ndarray, bound: float | None, *, below: bool):
    """Count values beyond one optional sanity bound."""
    if bound is None:
        return None
    comparison = values < bound if below else values > bound
    return int(np.count_nonzero(comparison))


def stage_band(output, index: int, definition: dict, inside: np.ndarray) -> dict:
    """Copy one source band exactly, applying only nodata normalization/masking."""
    source = open_raster(definition["path"])
    source_band = source.GetRasterBand(definition["source_band"])
    if source_band is None:
        raise ValueError(
            f"{relative(definition['path'])} has no band "
            f"{definition['source_band']}"
        )
    if source_band.GetDescription() != definition["source_description"]:
        raise ValueError(
            f"source description {source_band.GetDescription()!r} != "
            f"{definition['source_description']!r}"
        )
    source_dtype = gdal.GetDataTypeName(source_band.DataType)
    if source_dtype != definition["source_dtype"]:
        raise ValueError(
            f"{definition['name']} source dtype {source_dtype} != "
            f"{definition['source_dtype']}"
        )
    source_nodata = source_band.GetNoDataValue()
    if source_nodata is None or not math.isclose(
        source_nodata, definition["source_nodata"], rel_tol=0, abs_tol=0
    ):
        raise ValueError(
            f"{definition['name']} source nodata {source_nodata!r} != "
            f"{definition['source_nodata']!r}"
        )
    source_scale = source_band.GetScale()
    source_offset = source_band.GetOffset()
    if source_scale not in (None, 1.0) or source_offset not in (None, 0.0):
        raise ValueError(
            f"{definition['name']} has an unapplied scale/offset: "
            f"{source_scale!r}/{source_offset!r}"
        )

    values = source_band.ReadAsArray()
    if values is None or values.shape != (ROWS, COLS):
        raise RuntimeError(f"failed reading {definition['name']}")
    valid = valid_mask(values, source_nodata)
    valid_inside = valid & inside
    staged = np.full((ROWS, COLS), NODATA, dtype=np.float32)
    staged[valid_inside] = values[valid_inside]
    inside_values = values[valid_inside].astype(np.float64, copy=False)
    if inside_values.size == 0:
        raise ValueError(f"{definition['name']} has no valid cells inside Ethiopia")

    output_band = output.GetRasterBand(index)
    output_band.SetNoDataValue(NODATA)
    output_band.SetDescription(definition["name"])
    output_band.SetMetadataItem("unit", definition["unit"])
    output_band.SetMetadataItem("source", definition["source"])
    output_band.SetMetadataItem("source_file", relative(definition["path"]))
    output_band.SetMetadataItem("source_band", str(definition["source_band"]))
    output_band.SetMetadataItem("resample", "none")
    output_band.SetMetadataItem("logical_dtype", definition["logical_dtype"])
    output_band.SetMetadataItem(
        "overview_resample", definition["overview_resample"]
    )
    output_band.WriteArray(staged)

    statistics = {
        "valid_inside_cells": int(inside_values.size),
        "missing_inside_cells": int(np.count_nonzero(inside & ~valid)),
        "nonfinite_inside_cells": int(
            np.count_nonzero(inside & ~np.isfinite(values))
        ),
        "populated_outside_cells_before_mask": int(
            np.count_nonzero(~inside & valid)
        ),
        "min": float(inside_values.min()),
        "median": float(np.median(inside_values)),
        "max": float(inside_values.max()),
        "below_sanity_range_cells": range_count(
            inside_values, definition["valid_min"], below=True
        ),
        "above_sanity_range_cells": range_count(
            inside_values, definition["valid_max"], below=False
        ),
    }

    outside_values = values[~inside & valid].astype(np.float64, copy=False)
    statistics["outside_below_sanity_range_cells"] = range_count(
        outside_values, definition["valid_min"], below=True
    )
    statistics["outside_above_sanity_range_cells"] = range_count(
        outside_values, definition["valid_max"], below=False
    )

    if definition["logical_dtype"] == "uint8":
        rounded = np.rint(inside_values)
        noninteger_count = int(np.count_nonzero(inside_values != rounded))
        observed = sorted(int(value) for value in np.unique(rounded))
        illegal = sorted(set(observed) - set(LANDCOVER_CLASSES))
        statistics["noninteger_class_cells"] = noninteger_count
        statistics["observed_classes"] = observed
        statistics["class_counts"] = {
            str(code): int(np.count_nonzero(rounded == code))
            for code in observed
        }
        statistics["illegal_classes"] = illegal
        if noninteger_count or illegal:
            raise ValueError(
                f"landcover has noninteger_count={noninteger_count}, "
                f"illegal_classes={illegal}"
            )

    report = {
        "index": index,
        "name": definition["name"],
        "description": definition["description"],
        "group": definition["group"],
        "unit": definition["unit"],
        "physical_dtype": "float32",
        "logical_dtype": definition["logical_dtype"],
        "nodata": NODATA,
        "source": definition["source"],
        "source_file": relative(definition["path"]),
        "source_band": definition["source_band"],
        "source_band_description": definition["source_description"],
        "source_dtype": source_dtype,
        "source_nodata": source_nodata,
        "source_unit_tag": source_band.GetUnitType(),
        "source_scale": source_scale,
        "source_offset": source_offset,
        "resample": "none; exact-grid copy",
        "overview_resample": definition["overview_resample"],
        "valid_min": definition["valid_min"],
        "valid_max": definition["valid_max"],
        "notes": definition["notes"],
        "provenance": definition["provenance"],
        "statistics": statistics,
    }
    source = None
    return report


def validate_round_trip(
    definitions: list[dict], inside: np.ndarray, reports: list[dict]
) -> None:
    """Prove staged values and canonical descriptions survived exactly."""
    assert_matches_grid(OUTPUT_TEMP_PATH)
    output = open_raster(OUTPUT_TEMP_PATH)
    if output.RasterCount != len(definitions):
        raise ValueError(
            f"staged band count {output.RasterCount} != {len(definitions)}"
        )

    descriptions = []
    for index, (definition, report) in enumerate(
        zip(definitions, reports), start=1
    ):
        source = open_raster(definition["path"])
        source_band = source.GetRasterBand(definition["source_band"])
        source_values = source_band.ReadAsArray()
        valid = valid_mask(source_values, source_band.GetNoDataValue())
        expected = np.full((ROWS, COLS), NODATA, dtype=np.float32)
        expected[valid & inside] = source_values[valid & inside]

        output_band = output.GetRasterBand(index)
        output_values = output_band.ReadAsArray()
        mismatches = int(np.count_nonzero(output_values != expected))
        report["round_trip_mismatches"] = mismatches
        if mismatches:
            raise ValueError(
                f"{definition['name']} has {mismatches} source-copy mismatches"
            )
        descriptions.append(output_band.GetDescription())
        source = None

    expected_descriptions = [definition["name"] for definition in definitions]
    if descriptions != expected_descriptions:
        raise ValueError(
            f"staged descriptions {descriptions} != {expected_descriptions}"
        )
    output = None


def validate_rainfall_cv() -> dict:
    """Confirm each staged CV equals its supplied std divided by mean."""
    output = open_raster(OUTPUT_TEMP_PATH)
    descriptions = band_indices_by_name(output)
    checks = {}
    for season in SEASON_MONTHS:
        mean = output.GetRasterBand(
            descriptions[f"{season}_mean_mm"]
        ).ReadAsArray()
        standard_deviation = output.GetRasterBand(
            descriptions[f"{season}_std_mm"]
        ).ReadAsArray()
        coefficient = output.GetRasterBand(
            descriptions[f"{season}_cv"]
        ).ReadAsArray()
        valid = (
            np.isfinite(mean)
            & np.isfinite(standard_deviation)
            & np.isfinite(coefficient)
            & (mean != NODATA)
            & (standard_deviation != NODATA)
            & (coefficient != NODATA)
            & (mean > 0)
        )
        expected = standard_deviation[valid].astype(np.float64) / mean[
            valid
        ].astype(np.float64)
        observed = coefficient[valid].astype(np.float64)
        if not observed.size:
            raise ValueError(f"{season} CV has no valid positive-mean cells to compare")
        differences = np.abs(expected - observed)
        mismatches = int(
            np.count_nonzero(~np.isclose(expected, observed, rtol=1e-5, atol=1e-6))
        )
        if mismatches:
            raise ValueError(
                f"{season} CV has {mismatches} cells inconsistent with std/mean"
            )
        checks[season] = {
            "compared_cells": int(observed.size),
            "tolerance_mismatches": mismatches,
            "maximum_absolute_difference": float(differences.max()),
        }
    output = None
    return checks


def main() -> None:
    """Build and atomically publish the staged raster and its audit."""
    OUTPUT_TEMP_PATH.unlink(missing_ok=True)
    REPORT_TEMP_PATH.unlink(missing_ok=True)
    try:
        validate_source_layouts()
        fingerprints = source_fingerprints()
        rainfall_records = load_rainfall_metadata()
        definitions = build_band_definitions(rainfall_records)
        inside = read_mask()
        output = create_output(len(definitions))
        reports = [
            stage_band(output, index, definition, inside)
            for index, definition in enumerate(definitions, start=1)
        ]
        output.FlushCache()
        output = None

        validate_round_trip(definitions, inside, reports)
        rainfall_consistency = validate_rainfall_cv()
        if source_fingerprints() != fingerprints:
            raise RuntimeError("input files changed during the staging build")
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
            "band_order": [definition["name"] for definition in definitions],
            "bands": reports,
            "landcover_classes": {
                str(code): label for code, label in LANDCOVER_CLASSES.items()
            },
            "rainfall_companion": {
                "path": relative(RAINFALL_METADATA_PATH),
                "records_verified": len(rainfall_records),
                "monthly_files_read": 0,
            },
            "rainfall_cv_consistency": rainfall_consistency,
            "temperature_mapping": {
                definition["name"]: {
                    "source_band": definition["source_band"],
                    "source_description": definition["source_description"],
                }
                for definition in definitions if definition["group"] == "temperature"
            },
            "source_fingerprints": fingerprints,
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
