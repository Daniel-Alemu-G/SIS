#!/usr/bin/env python3
"""Assemble and accept the approved 38-band Phase 1 raster package.

Only prepared ingredients are consumed; nothing is downloaded or re-derived.
Phosphorus is AfricanSoilGrids in owner-confirmed mg/100kg and is copied as-is.
EthioSIS receives the approved no-interpolation window and country mask, never
an all-nutrients-present mask. Source gaps and diagnostic outliers are retained.

Run with the approved .venv from the repository root. This script lives in
pipeline/phase1_build; its working outputs stay beside the generated
data/interim artifacts. Only accepted, finished products are copied to
data/processed. --validate-only rechecks the delivered package without
regenerating anything.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr
from osgeo_utils.samples.validate_cloud_optimized_geotiff import validate as validate_cog
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config import grid, stack  # noqa: E402

HERE = ROOT / "data/interim/stack_v1"
PROCESSED = ROOT / "data/processed"
RAW_PATH = HERE / "common_stack_v1_source.tif"
COG_PATH = HERE / "common_stack_v1.tif"
MANIFEST_PATH = HERE / "band_manifest.json"
VALIDATION_PATH = HERE / "stack_v1_validation.json"
BUILD_REPORT_PATH = HERE / "stack_v1_build_report.json"
COG_TEMP = HERE / "common_stack_v1.tmp.tif"
LOOKUP_PATH = ROOT / "data/interim/admin_grid/woreda_id_lookup.csv"
ADMIN0_PATH = ROOT / "data/inputs/Admin boundaries/eth_admin0.geojson"
WOREDA_PATH = ROOT / "data/inputs/Admin boundaries/eth_admin3_woreda.geojson"
CAPITALS_PATH = ROOT / "data/inputs/Admin boundaries/eth_admincapitals.geojson"
DEM90_PATH = ROOT / "data/inputs/Elevation/ethiopia_copernicus_dem_90m.tif"
SLOPE90_PATH = ROOT / "data/inputs/Slope/slope_albers_90m.tif"
REPORT_PATHS = {
    "admin": ROOT / "data/interim/admin_grid/admin_grid_report.json",
    "ethiosis": ROOT / "data/interim/input_normalization/normalization_report.json",
    "slope": ROOT / "data/interim/slope_mean_deg/slope_mean_deg_report.json",
    "soil": ROOT / "data/interim/soil_elevation/soil_elevation_report.json",
    "climate": ROOT / "data/interim/landcover_climate/landcover_climate_report.json",
    "phosphorus": ROOT / "data/interim/phosphorus/phosphorus_report.json",
    "storage": ROOT / "data/interim/storage_prototype/prototype_report.json",
}
# This is the former location of this exact builder before the Phase 1 scripts
# moved into pipeline/. It is the only accepted legacy owner of stack v1.
LEGACY_BUILD_SCRIPTS = frozenset({"data/interim/stack_v1/build_stack_v1.py"})

gdal.UseExceptions()
ogr.UseExceptions()


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fingerprints(paths) -> dict:
    # Output mtimes/timestamps are deliberately absent: reruns must be stable.
    return {relative(path): {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted(set(paths))}


def input_snapshot() -> dict:
    return {relative(path): {"size_bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for path in sorted((ROOT / "data/inputs").rglob("*")) if path.is_file()}


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def validate_phosphorus_approval(report: dict) -> None:
    if report.get("status") != "passed" or report.get("stack_eligible") is not True:
        raise ValueError("phosphorus must have passed and have explicit owner inclusion approval")
    if report.get("publication_ready") is not True or report.get("publication_blockers"):
        raise ValueError("phosphorus still has a publication blocker")
    if report.get("source", {}).get("upstream_product") != "AfricanSoilGrids":
        raise ValueError("phosphorus source must be owner-confirmed AfricanSoilGrids")
    bands = report.get("bands", [])
    if len(bands) != 1 or bands[0].get("name") != "sg_p":
        raise ValueError("approved phosphorus report must describe only sg_p")
    if bands[0].get("unit") != "mg/100kg" or bands[0].get("transform_applied") != "none":
        raise ValueError("phosphorus must be used as-is in mg/100kg, without a transform")
    approval = report.get("owner_approval", {})
    if approval.get("product") != "AfricanSoilGrids" or approval.get("unit") != "mg/100kg":
        raise ValueError("phosphorus lacks the matching owner provenance authority")


def load_reports() -> dict:
    reports = {name: json.loads(path.read_text()) for name, path in REPORT_PATHS.items()}
    for name, report in reports.items():
        if report.get("status") != "passed":
            raise ValueError(f"prepared {name} audit has not passed")
        if report.get("sources_unchanged") is False or report.get("source_unchanged") is False:
            raise ValueError(f"prepared {name} audit reports changed source inputs")
        recorded = report.get("source_fingerprints")
        if recorded is not None and fingerprints(ROOT / path for path in recorded) != recorded:
            raise ValueError(f"prepared {name} sources no longer match their audited checksums")
    validate_phosphorus_approval(reports["phosphorus"])
    derivation = reports["slope"]["derivation"]
    if derivation["order"] != "convert each valid 90 m pixel, then average":
        raise ValueError("slope audit does not prove the approved conversion order")
    if reports["slope"]["unit"] != "degrees":
        raise ValueError("prepared slope is not in degrees")
    return reports


def build_definitions(reports: dict) -> list[dict]:
    definitions = {}

    def add(record: dict, staging_file: str, staging_name: str | None = None):
        record = dict(record)
        name = record["name"]
        if name in definitions:
            raise ValueError(f"duplicate prepared band: {name}")
        record["staging_file"] = staging_file
        record["staging_band_name"] = staging_name or name
        record.setdefault("source_band_name", record.get("source_band_description") or name)
        record.setdefault("transform_applied", "none")
        record.setdefault("valid_min", None)
        record.setdefault("valid_max", None)
        record.setdefault("notes", "")
        definitions[name] = record

    admin = reports["admin"]
    for name, source_path, description, group, minimum, maximum, missing in (
        ("ethiopia_mask", ADMIN0_PATH, "Ethiopia membership at the canonical cell centre", "mask", 0, 1,
         "No missing mask pixels: 0 outside Ethiopia, 1 inside; shared nodata tag is not used by mask pixels."),
        ("woreda_id", WOREDA_PATH, "Persisted woreda raster ID at the canonical cell centre", "admin", 0,
         max(int(row["raster_id"]) for row in read_lookup(LOOKUP_PATH).values()),
         "0 means no woreda inside Ethiopia; -9999 means outside Ethiopia. Positive IDs join the persisted lookup."),
    ):
        add({"name": name, "description": description, "group": group,
             "unit": stack.CORE_UNITS[name], "source": "Owner-approved Ethiopia administrative boundaries v04 (2025-01-01)",
             "source_file": relative(source_path), "source_band_name": None,
             "resample": admin["rasterization"], "valid_min": minimum, "valid_max": maximum,
             "missing_value_semantics": missing,
             "provenance": {"approval": admin["owner_approval"], "version": "v04", "valid_on": "2025-01-01"}},
            admin["mask" if name == "ethiopia_mask" else "woreda_id"]["path"])

    for group in ("soil", "climate", "phosphorus"):
        report = reports[group]
        for band in report["bands"]:
            add(band, report["output"])

    phosphorus = definitions["sg_p"]
    phosphorus.update({
        "valid_min": 70.0, "valid_max": 4000.0,
        "valid_range_role": "Provisional specification diagnostic ONLY; never a filter, clipping rule or scoring threshold.",
        "owner_approval": reports["phosphorus"]["owner_approval"],
        "provenance": reports["phosphorus"]["source"],
        "notes": "AfricanSoilGrids, not global SoilGrids. Preserve supplied mg/100kg values and gaps; no scaling or back-transformation.",
    })

    ethiosis = reports["ethiosis"]
    for source_name in (*stack.ETHIOSIS_NUTRIENTS, "confidence"):
        name = f"es_{source_name}"
        confidence = source_name == "confidence"
        add({"name": name, "description": "EthioSIS nearest-sample proximity confidence" if confidence else f"EthioSIS {source_name} supplied IDW estimate",
             "group": "ethiosis", "unit": stack.CORE_UNITS[name],
             "source": "EthioSIS supplied IDW estimates (nutrient units on source CSV/header basis)",
             "source_file": ethiosis["source"], "source_band_name": source_name,
             "resample": ethiosis["operation"], "valid_min": 0.0 if confidence else None,
             "valid_max": 1.0 if confidence else None,
             "missing_value_semantics": (
                 "0 is valid: no sample within the supplied 20 km support; confidence describes proximity, not all-nutrient availability. -9999 outside Ethiopia."
                 if confidence else
                 "-9999 for source gaps or outside Ethiopia; individual nutrients may be missing at positive confidence. No imputation or all-nutrient masking."
             ),
             "notes": "No extra transform of supplied IDW values. No es_n, es_p or es_oc is synthesized. Nutrient units are not embedded in the GeoTIFF.",
             "provenance": {"operation": ethiosis["operation"], "missingness_policy": ethiosis["missingness_policy"],
                            "source_shape_rows_cols": ethiosis["source_shape_rows_cols"],
                            "discarded_eastern_column": ethiosis["discarded_eastern_column"][source_name]}},
            ethiosis["output"], source_name)

    slope = reports["slope"]
    add({"name": "slope_mean_deg", "description": "Mean slope degrees: convert valid native 90 m percentages before common-grid averaging",
         "group": "terrain", "unit": "degrees", "source": "Copernicus DEM-derived 90 m Albers percentage slope",
         "source_file": slope["source_90m"]["path"], "source_band_name": "native_slope_percent",
         "resample": slope["derivation"]["resampling"], "transform_applied": slope["derivation"]["formula"],
         "valid_min": slope["output_statistics"]["sanity_range_deg"][0],
         "valid_max": slope["output_statistics"]["sanity_range_deg"][1],
         "notes": "The existing 1 km percentage raster is reference-only. Exact-click 90 m slope remains separate and must be CRS-transformed, sampled, then converted explicitly.",
         "provenance": {"derivation": slope["derivation"], "native_source": slope["source_90m"],
                        "reference_only": slope["existing_1km_percentage"]}}, slope["output"])

    if set(definitions) != set(stack.CORE_BAND_NAMES):
        raise ValueError("prepared ingredient names do not exactly match the 38-band core contract")
    ordered = []
    for index, name in enumerate(stack.CORE_BAND_NAMES, start=1):
        record = definitions[name]
        record.update({"index": index, "physical_dtype": "float32", "logical_dtype": stack.logical_dtype(name),
                       "nodata": grid.NODATA, "overview_resample": stack.overview_method(name),
                       "resolution_label": "Nominal 1 km canonical EPSG:4326 cell; not an exact 90 m click value."})
        record.setdefault("valid_range_role", "Specification sanity diagnostic only; report outliers without clipping. Null bounds mean no bound was supplied.")
        record.setdefault("missing_value_semantics", "-9999 for source missing values or outside Ethiopia; represent missing values as JSON null, never zero.")
        # Staging audits have their own indices and statistics; final indices
        # and statistics belong to the real assembled raster, not those audits.
        record.pop("statistics", None)
        ordered.append(record)
    return ordered


def read_lookup(path: Path) -> dict[int, dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    active = {}
    seen_pcodes = set()
    seen_ids = set()
    for row in rows:
        raster_id = int(row["raster_id"])
        pcode = row["adm3_pcode"]
        if raster_id <= 0 or raster_id > np.iinfo(np.uint16).max or raster_id in seen_ids or pcode in seen_pcodes:
            raise ValueError("persisted woreda lookup has invalid or duplicate identities")
        seen_ids.add(raster_id)
        seen_pcodes.add(pcode)
        if int(row["active"]):
            active[raster_id] = row
    if not active:
        raise ValueError("persisted woreda lookup has no active rows")
    return active


def validate_lookup_source(lookup: dict) -> dict:
    source = gdal.OpenEx(str(WOREDA_PATH), gdal.OF_VECTOR | gdal.OF_READONLY)
    rows = {row["adm3_pcode"]: row for row in lookup.values()}
    observed = set()
    for feature in source.GetLayer(0):
        pcode = feature.GetField("adm3_pcode")
        if pcode not in rows or pcode in observed:
            raise ValueError("active woreda lookup identities do not reconcile to source polygons")
        observed.add(pcode)
        row = rows[pcode]
        fields = {key: key for key in ("adm3_pcode", "adm3_name", "adm2_pcode", "adm2_name", "adm1_pcode", "adm1_name")}
        fields.update({"source_version": "version", "valid_on": "valid_on"})
        for lookup_field, source_field in fields.items():
            if row[lookup_field] != str(feature.GetField(source_field) or ""):
                raise ValueError(f"{pcode}: lookup {lookup_field} differs from the approved source")
    if observed != set(rows):
        raise ValueError("active woreda lookup contains identities absent from source polygons")
    source = None
    return {"matched_source_polygons": len(observed), "identity_and_parent_metadata": "passed"}


class RasterReaders:
    """Keep band lifetimes attached to read-only datasets; access by name."""

    def __init__(self):
        self.datasets = {}

    def open(self, path: Path):
        if path not in self.datasets:
            dataset = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
            if dataset is None:
                raise ValueError(f"cannot open raster read-only: {relative(path)}")
            self.datasets[path] = dataset
        return self.datasets[path]

    def band(self, path: Path, name: str, *, alias: bool = False):
        dataset = self.open(path)
        indices = stack.band_indices_by_name(dataset, unnamed_alias=name if alias else None)
        if name not in indices:
            raise ValueError(f"{relative(path)} has no band named {name!r}")
        return dataset.GetRasterBand(indices[name])

    def close(self):
        self.datasets.clear()


def masked_copy(values: np.ndarray, inside: np.ndarray, nodata: float | None) -> np.ndarray:
    if values.shape != inside.shape or inside.dtype != np.bool_:
        raise ValueError("source and boolean country mask must have the same shape")
    valid = np.isfinite(values)
    if nodata is not None and math.isfinite(nodata):
        valid &= values != nodata
    result = np.full(values.shape, grid.NODATA, dtype=np.float32)
    result[inside & valid] = values[inside & valid]
    return result


def bitwise_mismatches(actual: np.ndarray, expected: np.ndarray) -> int:
    if actual.shape != expected.shape or actual.dtype != np.float32 or expected.dtype != np.float32:
        raise ValueError("bitwise comparisons require matching float32 shapes")
    return int(np.count_nonzero(actual.view(np.uint32) != expected.view(np.uint32)))


def read_prepared(record: dict, readers: RasterReaders, inside: np.ndarray) -> np.ndarray:
    path = ROOT / record["staging_file"]
    dataset = readers.open(path)
    grid.assert_matches_grid(dataset)
    band = readers.band(path, record["staging_band_name"])
    if band.DataType != gdal.GDT_Float32:
        raise ValueError(f"{record['name']} staging is not physical float32")
    if record["name"] != "ethiopia_mask" and band.GetNoDataValue() != grid.NODATA:
        raise ValueError(f"{record['name']} staging has unexpected nodata")
    if record["name"] == "sg_p":
        if dataset.GetMetadataItem("stack_eligible") != "true" or dataset.GetMetadataItem("publication_ready") != "true":
            raise ValueError("staged phosphorus is not explicitly owner-approved")
    values = band.ReadAsArray()
    if record["name"] == "ethiopia_mask":
        if not np.all(np.isin(values, (0, 1))):
            raise ValueError("country mask contains missing/non-binary pixels")
        return values
    return masked_copy(values, inside, band.GetNoDataValue())


def original_expected(record: dict, readers: RasterReaders, inside: np.ndarray) -> np.ndarray:
    # Vector rasterization and slope conversion were already performed and
    # audited. Reuse those outputs; do not quietly rerun those derivations.
    if record["name"] in {"ethiopia_mask", "woreda_id", "slope_mean_deg"}:
        return read_prepared(record, readers, inside)
    path = ROOT / record["source_file"]
    dataset = readers.open(path)
    if record["group"] != "ethiosis":
        grid.assert_matches_grid(dataset)
    name = record["source_band_name"]
    band = readers.band(path, name, alias=dataset.RasterCount == 1 and not dataset.GetRasterBand(1).GetDescription())
    expected_dtype = gdal.GDT_Byte if record["name"] == "landcover" else gdal.GDT_Float32
    if band.DataType != expected_dtype:
        raise ValueError(f"{record['name']} source dtype differs from the inventoried product")
    nodata = band.GetNoDataValue()
    if record["group"] == "ethiosis":
        if nodata is None or not math.isnan(nodata):
            raise ValueError(f"{record['name']} source no longer uses inventoried NaN nodata")
    elif nodata != (0 if record["name"] == "landcover" else grid.NODATA):
        raise ValueError(f"{record['name']} source nodata differs from the inventoried product")
    if band.GetScale() not in (None, 1.0) or band.GetOffset() not in (None, 0.0):
        raise ValueError(f"{record['name']} has an unapplied source scale/offset")
    values = band.ReadAsArray(0, 0, grid.COLS, grid.ROWS)
    return masked_copy(values, inside, band.GetNoDataValue())


def verify_ethiosis_source(readers: RasterReaders, reports: dict) -> dict:
    path = ROOT / reports["ethiosis"]["source"]
    source = readers.open(path)
    if (source.RasterXSize, source.RasterYSize) != (grid.COLS + 1, grid.ROWS):
        raise ValueError("EthioSIS source no longer has the approved one-column mismatch")
    indices = stack.band_indices_by_name(source)
    expected_names = (*stack.ETHIOSIS_NUTRIENTS, "confidence")
    if tuple(indices) != expected_names:
        raise ValueError("EthioSIS source nutrient descriptions/order have changed")
    view = gdal.Translate("", source, format="VRT", srcWin=(0, 0, grid.COLS, grid.ROWS))
    grid.assert_matches_grid(view)
    view = None
    details = {}
    for name in expected_names:
        band = readers.band(path, name)
        if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() is None or not math.isnan(band.GetNoDataValue()):
            raise ValueError(f"EthioSIS {name} no longer has the inventoried float32/NaN encoding")
        edge = band.ReadAsArray(grid.COLS, 0, 1, grid.ROWS)
        if name == "confidence":
            if not np.all(np.isfinite(edge) & (edge == 0)):
                raise ValueError("discarded EthioSIS confidence column is not finite zero")
        elif np.any(np.isfinite(edge) & (edge != grid.NODATA)):
            raise ValueError(f"discarded EthioSIS column has valid {name} values")
        details[name] = {"finite_cells": int(np.isfinite(edge).sum()), "nonzero_finite_cells": int(np.count_nonzero(np.isfinite(edge) & (edge != 0)))}
    confidence = readers.band(path, "confidence").ReadAsArray(0, 0, grid.COLS, grid.ROWS)
    if not np.all(np.isfinite(confidence) & (confidence >= 0) & (confidence <= 1)):
        raise ValueError("the approved EthioSIS confidence window is not fully valid in 0..1")
    return details


def band_statistics(values: np.ndarray, inside: np.ndarray, record: dict) -> dict:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{record['name']} contains non-finite physical output values")
    valid = values != grid.NODATA
    selected = values[valid] if record["name"] == "ethiopia_mask" else values[inside & valid]
    if not selected.size:
        raise ValueError(f"{record['name']} has no valid data")
    lower, upper = record["valid_min"], record["valid_max"]
    return {"statistics_domain": "whole grid" if record["name"] == "ethiopia_mask" else "inside Ethiopia",
            "valid_inside_cells": int(np.count_nonzero(inside & valid)),
            "missing_inside_cells": int(np.count_nonzero(inside & ~valid)),
            "nonfinite_cells": 0, "min": float(selected.min()), "median": float(np.median(selected)),
            "mean": float(selected.mean(dtype=np.float64)), "max": float(selected.max()),
            "below_sanity_range_cells": None if lower is None else int(np.count_nonzero(selected < lower)),
            "above_sanity_range_cells": None if upper is None else int(np.count_nonzero(selected > upper))}


def validate_logical_values(name: str, values: np.ndarray, lookup_ids: set[int], landcover_codes: set[int]) -> None:
    valid = values[values != grid.NODATA]
    if name == "ethiopia_mask":
        if np.any(values == grid.NODATA) or not np.all(np.isin(values, (0, 1))):
            raise ValueError("mask must have only meaningful 0/1 pixels")
    elif name in {"woreda_id", "landcover"}:
        if not np.all(valid == np.rint(valid)):
            raise ValueError(f"{name} has fractional categorical/ID values")
        allowed = lookup_ids | {0} if name == "woreda_id" else landcover_codes
        if not set(int(value) for value in np.unique(valid)) <= allowed:
            raise ValueError(f"{name} contains illegal categorical/ID values")


def configure_output(dataset) -> None:
    dataset.SetGeoTransform(grid.TRANSFORM.to_gdal())
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(grid.CRS)
    dataset.SetSpatialRef(spatial_ref)
    dataset.SetMetadataItem("AREA_OR_POINT", "Area")
    dataset.SetMetadataItem("stack_version", str(stack.STACK_VERSION))
    dataset.SetMetadataItem("role", "Approved Phase 1 raster package; scoring/access bands are absent")


def apply_overview_policy(dataset, indices: dict[str, int], factors: list[int]) -> None:
    if dataset.BuildOverviews("AVERAGE", factors) != gdal.CE_None:
        raise ValueError("could not create common average overview levels")
    # GTiff may cascade average overviews. Around nodata, an average of
    # averages loses contributor weights; non-divisible dimensions also change
    # support. Regenerate EVERY level directly from each full-resolution band.
    # This affects only display overviews, never original/full-resolution data.
    for name in indices:
        band = dataset.GetRasterBand(indices[name])
        method = stack.overview_method(name).upper()
        for index in range(band.GetOverviewCount()):
            if gdal.RegenerateOverview(band, band.GetOverview(index), method) != gdal.CE_None:
                raise ValueError(f"could not regenerate {name} {method.lower()} overview")


def assemble(definitions: list[dict], inside: np.ndarray, readers: RasterReaders) -> None:
    output = gdal.GetDriverByName("GTiff").Create(str(RAW_PATH), grid.COLS, grid.ROWS, len(definitions), gdal.GDT_Float32,
        options=["TILED=YES", f"BLOCKXSIZE={stack.BLOCK_SIZE}", f"BLOCKYSIZE={stack.BLOCK_SIZE}",
                 "COMPRESS=DEFLATE", "PREDICTOR=3", "INTERLEAVE=BAND", "NUM_THREADS=2"])
    configure_output(output)
    lookup_ids = set(read_lookup(LOOKUP_PATH))
    classes = set(int(value) for value in load_reports()["climate"]["landcover_classes"])
    for record in definitions:
        name = record["name"]
        values = read_prepared(record, readers, inside)
        mismatches = bitwise_mismatches(values, original_expected(record, readers, inside))
        if mismatches:
            raise ValueError(f"{name} has {mismatches} prepared/original copy mismatches")
        validate_logical_values(name, values, lookup_ids, classes)
        record["statistics"] = band_statistics(values, inside, record)
        record["original_source_copy_mismatches"] = mismatches
        band = output.GetRasterBand(record["index"])
        band.SetDescription(name)
        band.SetNoDataValue(grid.NODATA)
        band.SetUnitType(record["unit"])
        for key in ("unit", "source", "source_file", "logical_dtype", "physical_dtype", "overview_resample",
                    "resample", "transform_applied", "missing_value_semantics", "resolution_label", "staging_file", "staging_band_name"):
            band.SetMetadataItem(key, str(record[key]))
        band.WriteArray(values)
        print(f"Copied {name}: {record['statistics']['valid_inside_cells']:,} valid inside cells; no value transform during assembly.", flush=True)
    output.FlushCache()
    apply_overview_policy(output, stack.band_indices_by_name(output), stack.overview_factors(grid.COLS, grid.ROWS))
    output.FlushCache()
    output = None
    cog = gdal.Translate(str(COG_TEMP), str(RAW_PATH), format="COG", creationOptions=[
        "COMPRESS=DEFLATE", "PREDICTOR=FLOATING_POINT", "OVERVIEW_COMPRESS=DEFLATE",
        "OVERVIEW_PREDICTOR=FLOATING_POINT", f"BLOCKSIZE={stack.BLOCK_SIZE}",
        "OVERVIEWS=FORCE_USE_EXISTING", "NUM_THREADS=2"])
    if cog is None:
        raise ValueError("COG translation failed")
    cog = None
    COG_TEMP.replace(COG_PATH)


def verify_overviews(output_band, raw_band, record: dict, factors: list[int], lookup_ids: set[int], classes: set[int]) -> list[dict]:
    # Independently generate a single-band reference from full-resolution
    # values. Compare its actual pixels, not merely an overview-resample tag.
    reference = gdal.GetDriverByName("MEM").Create("", grid.COLS, grid.ROWS, 1, gdal.GDT_Float32)
    source = reference.GetRasterBand(1)
    source.SetNoDataValue(grid.NODATA)
    source.WriteArray(output_band.ReadAsArray())
    details = []
    for index, factor in enumerate(factors):
        actual_overview = output_band.GetOverview(index)
        target = gdal.GetDriverByName("MEM").Create("", actual_overview.XSize, actual_overview.YSize, 1, gdal.GDT_Float32)
        target_band = target.GetRasterBand(1)
        target_band.SetNoDataValue(grid.NODATA)
        if gdal.RegenerateOverview(source, target_band, record["overview_resample"].upper()) != gdal.CE_None:
            raise ValueError(f"{record['name']}: could not generate independent overview reference")
        actual = actual_overview.ReadAsArray()
        expected = target_band.ReadAsArray()
        method_mismatches = bitwise_mismatches(actual, expected)
        copied_mismatches = bitwise_mismatches(actual, raw_band.GetOverview(index).ReadAsArray())
        if method_mismatches or copied_mismatches:
            raise ValueError(f"{record['name']} overview {factor}: method={method_mismatches}, copied={copied_mismatches} mismatches")
        if not np.all(np.isfinite(actual)):
            raise ValueError(f"{record['name']} overview contains non-finite physical values")
        validate_logical_values(record["name"], actual, lookup_ids, classes)
        details.append({"factor": factor, "shape_rows_cols": list(actual.shape), "method": record["overview_resample"],
                        "reference_mismatches": method_mismatches, "cog_copy_mismatches": copied_mismatches})
        target = None
    reference = None
    return details


def json_value(value) -> float | None:
    return None if not np.isfinite(value) or value == grid.NODATA else float(value)


def sample_native(path: Path, alias: str, lon: float, lat: float, readers: RasterReaders) -> dict:
    dataset = readers.open(path)
    transform = Transformer.from_crs(grid.CRS, dataset.GetSpatialRef().ExportToWkt(), always_xy=True)
    x, y = transform.transform(lon, lat)
    col_float, row_float = gdal.ApplyGeoTransform(gdal.InvGeoTransform(dataset.GetGeoTransform()), x, y)
    row, col = math.floor(row_float), math.floor(col_float)
    if not (0 <= row < dataset.RasterYSize and 0 <= col < dataset.RasterXSize):
        raise ValueError(f"spot check is outside {relative(path)} after CRS conversion")
    band = readers.band(path, alias, alias=True)
    value = float(band.ReadAsArray(col, row, 1, 1)[0, 0])
    if not np.isfinite(value) or value == band.GetNoDataValue():
        raise ValueError(f"spot check has no valid value in {relative(path)}")
    return {"file": relative(path), "row": row, "col": col, "native_x": x, "native_y": y,
            "value": value, "resolution_label": "Exact source pixel, nominal 90 m; not the 1 km mean."}


def original_point_value(record: dict, readers: RasterReaders, row: int, col: int) -> float | None:
    prepared = record["name"] in {"ethiopia_mask", "woreda_id", "slope_mean_deg"}
    path = ROOT / record["staging_file" if prepared else "source_file"]
    name = record["staging_band_name" if prepared else "source_band_name"]
    dataset = readers.open(path)
    band = readers.band(path, name, alias=dataset.RasterCount == 1 and not dataset.GetRasterBand(1).GetDescription())
    value = band.ReadAsArray(col, row, 1, 1)[0, 0]
    return None if not np.isfinite(value) or value == band.GetNoDataValue() else float(value)


def geographic_checks(dataset, indices: dict, records: dict, readers: RasterReaders, inside: np.ndarray, lookup: dict) -> dict:
    capitals = json.loads(CAPITALS_PATH.read_text())["features"]
    locations = {}
    for name in ("Addis Ababa", "Bahir Dar", "Godey"):
        matches = [feature for feature in capitals if feature["properties"]["name"] == name]
        if len(matches) != 1:
            raise ValueError(f"expected one supplied administrative capital for diagnostic site {name}")
        locations[name] = matches[0]["geometry"]["coordinates"][:2]
    # A labelled representative depression site, not a new named-object input.
    locations["Danakil/Dallol diagnostic site"] = [40.3, 14.24]
    results = {}
    source_woredas = gdal.OpenEx(str(WOREDA_PATH), gdal.OF_VECTOR | gdal.OF_READONLY)
    geometries = {feature.GetField("adm3_pcode"): feature.GetGeometryRef().Clone() for feature in source_woredas.GetLayer(0)}
    for name, (lon, lat) in locations.items():
        col_float, row_float = gdal.ApplyGeoTransform(gdal.InvGeoTransform(grid.TRANSFORM.to_gdal()), lon, lat)
        row, col = math.floor(row_float), math.floor(col_float)
        if not (0 <= row < grid.ROWS and 0 <= col < grid.COLS) or not inside[row, col]:
            raise ValueError(f"geographic diagnostic site {name} is outside the approved Ethiopia mask")
        values = {band_name: json_value(dataset.GetRasterBand(index).ReadAsArray(col, row, 1, 1)[0, 0]) for band_name, index in indices.items()}
        centre = gdal.ApplyGeoTransform(grid.TRANSFORM.to_gdal(), col + 0.5, row + 0.5)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint_2D(*centre)
        woreda = lookup[int(values["woreda_id"])]
        if not geometries[woreda["adm3_pcode"]].Intersects(point):
            raise ValueError(f"{name}: persisted raster ID does not identify the polygon containing the cell centre")
        source_comparisons = {}
        for band_name, record in records.items():
            source_value = original_point_value(record, readers, row, col)
            if source_value != values[band_name]:
                raise ValueError(f"{name}: {band_name} source/stack spot mismatch")
            source_comparisons[band_name] = {"stack": values[band_name], "source_or_approved_derived": source_value}
        elevation90 = sample_native(DEM90_PATH, "native_elevation_m", lon, lat, readers)
        slope90 = sample_native(SLOPE90_PATH, "native_slope_percent", lon, lat, readers)
        slope90["native_unit"] = "percent"
        slope90["degrees_after_sampling"] = float(np.degrees(np.arctan(slope90["value"] / 100)))
        results[name] = {"lon": lon, "lat": lat, "coordinate_provenance": "Supplied administrative capital" if name != "Danakil/Dallol diagnostic site" else "Representative depression test coordinate only",
                         "common_grid_row": row, "common_grid_col": col, "cell_centre_lon_lat": list(centre),
                         "resolution_label": "1 km cell values; separate exact native-pixel values below",
                         "woreda_cell_context": {key: woreda[key] for key in ("adm3_pcode", "adm3_name", "adm2_name", "adm1_name")},
                         "cell_centre_polygon_check": "passed", "source_round_trips": source_comparisons,
                         "elevation90_m": elevation90, "slope90_percent_and_converted_degrees": slope90}
    checks = {"danakil_below_sea_level": results["Danakil/Dallol diagnostic site"]["source_round_trips"]["elevation_m"]["stack"] < 0,
              "highland_rain_exceeds_somali_rain": results["Bahir Dar"]["source_round_trips"]["kiremt_mean_mm"]["stack"] > results["Godey"]["source_round_trips"]["kiremt_mean_mm"]["stack"]}
    if not all(checks.values()):
        raise ValueError(f"geographic acceptance failed: {checks}")
    source_woredas = None
    return {"checks": checks, "addis_rough_elevation_target_m": 2350, "addis_target_role": "Rough specification diagnostic, not a clipping limit; report the actual sampled elevation.", "sites": results}


def validate_package(path: Path, manifest: dict, *, lookup_path: Path = LOOKUP_PATH) -> dict:
    before = input_snapshot()
    reports = load_reports()
    current_inputs = fingerprints(ROOT / name for name in manifest["build_inputs"])
    if current_inputs != manifest["build_inputs"]:
        raise ValueError("source, prepared ingredient, audit or build-code checksum changed since assembly")
    if sha256(path) != manifest["raster_sha256"]:
        raise ValueError("COG checksum does not match manifest")
    readers = RasterReaders()
    output = readers.open(path)
    indices = stack.assert_stack_matches_manifest(output, manifest)
    warnings, errors, _ = validate_cog(output, full_check=True)
    if errors:
        raise ValueError(f"full GDAL COG structural validation failed: {errors}")
    raw = readers.open(RAW_PATH)
    raw_indices = stack.band_indices_by_name(raw)
    grid.assert_matches_grid(raw)
    records = stack.validate_manifest(manifest)
    mask = output.GetRasterBand(indices["ethiopia_mask"]).ReadAsArray()
    inside = mask == 1
    validate_logical_values("ethiopia_mask", mask, set(), set())
    lookup = read_lookup(lookup_path)
    if sha256(lookup_path) != manifest["lookup"]["sha256"]:
        raise ValueError("woreda lookup checksum does not match the manifest")
    lookup_source_check = validate_lookup_source(lookup)
    classes = set(int(value) for value in reports["climate"]["landcover_classes"])
    confidence = output.GetRasterBand(indices["es_confidence"]).ReadAsArray()
    if np.any(inside & ((confidence == grid.NODATA) | ~np.isfinite(confidence) | (confidence < 0) | (confidence > 1))):
        raise ValueError("confidence must be valid in 0..1 in every inside-country cell")
    discarded = verify_ethiosis_source(readers, reports)
    bands = []
    for name, record in records.items():
        band = output.GetRasterBand(indices[name])
        values = band.ReadAsArray()
        if name != "ethiopia_mask" and np.any(values[~inside] != grid.NODATA):
            raise ValueError(f"{name} has non-nodata pixels outside Ethiopia")
        validate_logical_values(name, values, set(lookup), classes)
        original_mismatches = bitwise_mismatches(values, original_expected(record, readers, inside))
        prepared_mismatches = bitwise_mismatches(values, read_prepared(record, readers, inside))
        if original_mismatches or prepared_mismatches:
            raise ValueError(f"{name}: final source/prepared bitwise comparison failed")
        statistics = band_statistics(values, inside, record)
        if statistics != record["statistics"]:
            raise ValueError(f"{name}: manifest statistics differ from actual raster")
        details = {"name": name, "statistics": statistics, "original_source_mismatches": original_mismatches,
                   "prepared_copy_mismatches": prepared_mismatches, "outside_country_check": "passed",
                   "overviews": verify_overviews(band, raw.GetRasterBand(raw_indices[name]), record, manifest["overview_factors"], set(lookup), classes)}
        if name.startswith("es_") and name != "es_confidence":
            valid = values != grid.NODATA
            invalid_at_zero = int(np.count_nonzero(inside & valid & (confidence == 0)))
            if invalid_at_zero:
                raise ValueError(f"{name} has valid nutrient values at zero confidence")
            details["valid_with_zero_confidence"] = invalid_at_zero
            details["missing_with_positive_confidence"] = int(np.count_nonzero(inside & ~valid & (confidence > 0)))
        bands.append(details)
    rainfall = {}
    for season in ("kiremt", "belg", "bega"):
        mean, std, cv = [output.GetRasterBand(indices[f"{season}_{suffix}"]).ReadAsArray()[inside].astype(np.float64)
                         for suffix in ("mean_mm", "std_mm", "cv")]
        if np.any(mean <= 0) or np.any(std < 0) or np.any(cv < 0):
            raise ValueError(f"{season}: invalid supplied rainfall mean/std/CV")
        # Numerical tolerances for float32 storage, not scientific thresholds.
        error = np.abs(cv - std / mean)
        mismatches = int(np.count_nonzero(~np.isclose(cv, std / mean, rtol=1e-5, atol=1e-6)))
        if mismatches:
            raise ValueError(f"{season}: CV does not match supplied std/mean")
        rainfall[season] = {"compared_cells": int(mean.size), "mismatches": mismatches, "max_absolute_error": float(error.max())}
    geography = geographic_checks(output, indices, records, readers, inside, lookup)
    woreda_values = output.GetRasterBand(indices["woreda_id"]).ReadAsArray()
    observed_ids = set(int(value) for value in np.unique(woreda_values[inside]) if value > 0)
    zero_cell_ids = sorted(set(lookup) - observed_ids)
    if zero_cell_ids != sorted(row["raster_id"] for row in reports["admin"]["woreda_id"]["zero_cell_woredas"]):
        raise ValueError("zero-cell woreda findings no longer match the approved rasterization audit")
    readers.close()
    if input_snapshot() != before:
        raise ValueError("an input file's inventory, size or modification time changed during validation")
    return {"status": "passed", "scope": "Phase 1 raster acceptance only; not complete Phase 1/PostGIS/PMTiles sign-off",
            "raster": manifest["raster_file"], "raster_sha256": manifest["raster_sha256"], "stack_version": stack.STACK_VERSION,
            "core_bands": len(records), "grid_check": "passed", "cog_full_structural_check": "passed", "cog_warnings": warnings,
            "inside_country_cells": int(inside.sum()), "outside_country_cells": int(inside.size - inside.sum()),
            "all_original_inputs_unchanged": True, "original_input_files": len(before),
            "bands": bands, "rainfall_cv_consistency": rainfall, "discarded_ethiosis_eastern_column": discarded,
            "woreda_lookup": {"rows": len(lookup), "observed_positive_ids": len(observed_ids), "zero_cell_ids": zero_cell_ids,
                              "source": manifest["lookup"]["file"], "sha256": manifest["lookup"]["sha256"],
                              "source_reconciliation": lookup_source_check},
            "geographic_checks": geography,
            "skipped_sanity_bounds": [{"name": name, "lower_bound": record["valid_min"], "upper_bound": record["valid_max"],
                                       "reason": "No corresponding scientific/sanity bound supplied; no value invented."}
                                      for name, record in records.items() if record["valid_min"] is None or record["valid_max"] is None],
            "remaining_phase1": ["Disposable PostGIS target/tooling and named-object loaders/index/query acceptance", "roads.pmtiles and HTTP Range/display smoke test", "Complete Phase 1 sign-off"],
            "deferred_optional_bands": ["dist_main_road_km", "dist_any_road_km", "dist_town_km"]}


def build_manifest(definitions: list[dict], build_inputs: dict, reports: dict) -> dict:
    manifest = {"stack_version": stack.STACK_VERSION, "grid": stack.grid_metadata(),
                "raster_file": relative(PROCESSED / COG_PATH.name), "raster_sha256": sha256(COG_PATH),
                "raster_size_bytes": COG_PATH.stat().st_size,
                "overview_factors": stack.overview_factors(grid.COLS, grid.ROWS), "bands": definitions,
                "storage": {"physical_dtype": "float32", "nodata": grid.NODATA, "compression": "DEFLATE", "predictor": 3,
                            "block_size": [stack.BLOCK_SIZE, stack.BLOCK_SIZE], "continuous_overviews": "average",
                            "categorical_id_mask_overviews": "mode",
                            "overview_generation": "Regenerate every level directly from full-resolution values, not an average-of-averages pyramid.",
                            "overview_edge_policy": "Average valid contributors; the full-resolution Ethiopia mask/nodata policy is authoritative."},
                "build_script": relative(Path(__file__).resolve()), "build_inputs": build_inputs,
                "runtime": {"python": platform.python_version(), "gdal": gdal.VersionInfo("RELEASE_NAME"), "numpy": np.__version__},
                "lookup": {"file": relative(PROCESSED / LOOKUP_PATH.name), "sha256": sha256(LOOKUP_PATH), "rows": len(read_lookup(LOOKUP_PATH))},
                "assumptions": ["All grid values are imported from pipeline/grid.yaml; no input is silently regridded.",
                                "AfricanSoilGrids phosphorus is owner-approved as-is in mg/100kg; no scaling, back-transform, clipping or imputation.",
                                "sg_n is the sole nitrogen input, in mg/kg without rescaling. No EthioSIS N/P/OC raster inputs exist.",
                                "Source gaps are retained, including nutrient-specific EthioSIS gaps at positive confidence.",
                                "1 km slope averages degrees computed at valid native 90 m pixels. The existing 1 km percent raster is not a substitute.",
                                "Rainfall reuses the supplied seasonal raster, means/std/CV and companion season-years; native CHIRPS is about 5 km, not new 1 km detail.",
                                "Land cover is a dominant class, not cropland percentage; categorical/ID/mask overviews never average codes.",
                                "Temperature is the supplied product in degrees Celsius with explicit source-description mappings; no unreported producer/version is invented.",
                                "90 m elevation and Albers percentage-slope inputs remain separate unchanged rasters. Exact slope sampling must transform CRS before explicit conversion.",
                                "Mask/admin bands use the approved v04 polygons and cell-centre rasterization. Two small woredas have no canonical cell centre.",
                                "Sanity ranges are diagnostics, not scoring thresholds; outliers are preserved. No scoring/API/UI or named-object infrastructure is implemented."],
                "owner_phosphorus_approval": reports["phosphorus"]["owner_approval"],
                "deferred_optional_bands": ["dist_main_road_km", "dist_any_road_km", "dist_town_km"],
                "phase1_scope": "Finished core raster package; PostGIS and PMTiles delivery is still owner-deferred."}
    stack.validate_manifest(manifest)
    return manifest


def publish(manifest: dict) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    final_cog = PROCESSED / COG_PATH.name
    final_manifest = PROCESSED / MANIFEST_PATH.name
    previous = None
    if final_cog.exists() or final_manifest.exists():
        if not final_cog.is_file() or not final_manifest.is_file():
            raise ValueError("existing final stack/manifest pair is incomplete; do not overwrite unfamiliar work")
        previous = json.loads(final_manifest.read_text())
        known_builders = {manifest["build_script"], *LEGACY_BUILD_SCRIPTS}
        if previous.get("build_script") not in known_builders:
            raise ValueError("existing final stack was not built by this step; do not overwrite it")
        if sha256(final_cog) != manifest["raster_sha256"]:
            raise ValueError("existing stack v1 has different content; a published raster change requires an explicit version decision")
    sources = (COG_PATH, LOOKUP_PATH, VALIDATION_PATH, MANIFEST_PATH)
    # Check ownership of every target before replacing any finished product.
    for source in sources:
        destination = PROCESSED / source.name
        if destination.exists() and source not in (COG_PATH, MANIFEST_PATH):
            known_builders = {manifest["build_script"], *LEGACY_BUILD_SCRIPTS}
            if previous is None or previous.get("build_script") not in known_builders:
                raise ValueError(f"refusing to replace unfamiliar deliverable {relative(destination)}")
    # A file is copied only after acceptance. Publish the matching manifest last.
    for source in sources:
        destination = PROCESSED / source.name
        temporary = destination.with_name(destination.name + ".tmp")
        shutil.copyfile(source, temporary)
        if sha256(temporary) != sha256(source):
            raise ValueError(f"finished-product copy checksum mismatch: {relative(destination)}")
        temporary.replace(destination)
    stack.assert_stack_matches_manifest(final_cog, json.loads(final_manifest.read_text()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true", help="Read-only acceptance of the delivered package, with no regeneration")
    args = parser.parse_args()
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetCacheMax(768 * 1024 * 1024)
    if args.validate_only:
        manifest = json.loads((PROCESSED / MANIFEST_PATH.name).read_text())
        result = validate_package(PROCESSED / COG_PATH.name, manifest, lookup_path=PROCESSED / LOOKUP_PATH.name)
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        return
    HERE.mkdir(parents=True, exist_ok=True)
    before = input_snapshot()
    reports = load_reports()
    definitions = build_definitions(reports)
    readers = RasterReaders()
    mask_record = definitions[0]
    mask_dataset = readers.open(ROOT / mask_record["staging_file"])
    grid.assert_matches_grid(mask_dataset)
    mask = readers.band(ROOT / mask_record["staging_file"], "ethiopia_mask").ReadAsArray()
    validate_logical_values("ethiopia_mask", mask, set(), set())
    inside = mask == 1
    verify_ethiosis_source(readers, reports)
    input_paths = set(REPORT_PATHS.values()) | {LOOKUP_PATH, CAPITALS_PATH, DEM90_PATH, SLOPE90_PATH,
        ROOT / "data/inputs/Rainfall/seasonal/rainfall_seasonal_1km.json", ROOT / "pipeline/grid.yaml",
        ROOT / "config/grid.py", ROOT / "config/stack.py", Path(__file__).resolve()}
    for record in definitions:
        input_paths.update({ROOT / record["source_file"], ROOT / record["staging_file"]})
    for report in reports.values():
        input_paths.update(ROOT / name for name in report.get("source_fingerprints", {}))
    for step in ("build_admin_grid.py", "normalize_ethiosis.py", "build_slope_mean_deg.py",
                 "stage_soil_elevation.py", "stage_landcover_climate.py", "stage_phosphorus.py"):
        input_paths.add(ROOT / "pipeline/phase1_build" / step)
    original_fingerprints = fingerprints(input_paths)
    assemble(definitions, inside, readers)
    readers.close()
    manifest = build_manifest(definitions, original_fingerprints, reports)
    write_json(MANIFEST_PATH, manifest)
    validation = validate_package(COG_PATH, manifest)
    if input_snapshot() != before or fingerprints(input_paths) != original_fingerprints:
        raise ValueError("original or prepared inputs changed during assembly/acceptance")
    write_json(VALIDATION_PATH, validation)
    publish(manifest)
    build_report = {"status": "passed", "stack_version": stack.STACK_VERSION, "core_bands": len(definitions),
                    "all_original_inputs_unchanged": True, "original_input_files": len(before),
                    "outputs": fingerprints([PROCESSED / source.name for source in (COG_PATH, MANIFEST_PATH, LOOKUP_PATH, VALIDATION_PATH)]),
                    "assembly_reference": relative(RAW_PATH), "validation": relative(PROCESSED / VALIDATION_PATH.name),
                    "scope": manifest["phase1_scope"]}
    write_json(BUILD_REPORT_PATH, build_report)
    print(json.dumps(build_report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
