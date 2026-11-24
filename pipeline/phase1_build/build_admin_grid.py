#!/usr/bin/env python3
"""Build the approved Ethiopia mask and deterministic woreda-ID grid."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr

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

ADMIN0_PATH = ROOT / "data/inputs/Admin boundaries/eth_admin0.geojson"
WOREDA_PATH = ROOT / "data/inputs/Admin boundaries/eth_admin3_woreda.geojson"
HERE = ROOT / "data/interim/admin_grid"
MASK_PATH = HERE / "ethiopia_mask.tif"
WOREDA_ID_PATH = HERE / "woreda_id.tif"
LOOKUP_PATH = HERE / "woreda_id_lookup.csv"
REPORT_PATH = HERE / "admin_grid_report.json"
MASK_TEMP = HERE / "ethiopia_mask.tmp.tif"
WOREDA_TEMP = HERE / "woreda_id.tmp.tif"
LOOKUP_TEMP = HERE / "woreda_id_lookup.tmp.csv"

LOOKUP_FIELDS = (
    "raster_id",
    "adm3_pcode",
    "adm3_name",
    "adm2_pcode",
    "adm2_name",
    "adm1_pcode",
    "adm1_name",
    "source_version",
    "valid_on",
    "active",
)

gdal.UseExceptions()
ogr.UseExceptions()
gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")


def remove_own_output(path: Path) -> None:
    if path.exists():
        path.unlink()


def source_snapshot(path: Path) -> dict:
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def expected_spatial_ref():
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(CRS)
    spatial_ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return spatial_ref


def open_layer(path: Path):
    dataset = gdal.OpenEx(str(path), gdal.OF_VECTOR | gdal.OF_READONLY)
    if dataset is None or dataset.GetLayerCount() != 1:
        raise RuntimeError(f"expected one readable vector layer in {path}")
    layer = dataset.GetLayer(0)
    actual_crs = layer.GetSpatialRef()
    if actual_crs is None:
        raise RuntimeError(f"{path} has no CRS")
    actual_crs = actual_crs.Clone()
    actual_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if not actual_crs.IsSame(expected_spatial_ref()):
        raise RuntimeError(f"{path} is not equivalent to {CRS}")
    return dataset, layer


def read_woredas(layer) -> list[dict]:
    required = (
        "adm3_pcode",
        "adm3_name",
        "adm2_pcode",
        "adm2_name",
        "adm1_pcode",
        "adm1_name",
        "version",
        "valid_on",
    )
    available = {
        layer.GetLayerDefn().GetFieldDefn(index).GetName()
        for index in range(layer.GetLayerDefn().GetFieldCount())
    }
    missing = sorted(set(required) - available)
    if missing:
        raise RuntimeError(f"woreda source lacks fields: {missing}")

    records = []
    seen_pcodes = set()
    layer.ResetReading()
    for feature in layer:
        pcode = feature.GetField("adm3_pcode")
        geometry = feature.GetGeometryRef()
        if not pcode or pcode in seen_pcodes:
            raise RuntimeError(f"missing or duplicate adm3_pcode: {pcode!r}")
        if geometry is None or geometry.IsEmpty() or not geometry.IsValid():
            raise RuntimeError(f"invalid woreda geometry for {pcode}")
        seen_pcodes.add(pcode)
        records.append(
            {
                "adm3_pcode": pcode,
                "adm3_name": feature.GetField("adm3_name") or "",
                "adm2_pcode": feature.GetField("adm2_pcode") or "",
                "adm2_name": feature.GetField("adm2_name") or "",
                "adm1_pcode": feature.GetField("adm1_pcode") or "",
                "adm1_name": feature.GetField("adm1_name") or "",
                "source_version": feature.GetField("version") or "",
                "valid_on": str(feature.GetField("valid_on") or ""),
                "geometry": geometry.Clone(),
            }
        )
    return records


def assign_ids(records: list[dict]) -> list[dict]:
    existing = {}
    historical_rows = {}
    if LOOKUP_PATH.exists():
        with LOOKUP_PATH.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != LOOKUP_FIELDS:
                raise RuntimeError("existing lookup schema is not recognized")
            for row in reader:
                pcode = row["adm3_pcode"]
                raster_id = int(row["raster_id"])
                if pcode in existing or raster_id in existing.values():
                    raise RuntimeError("existing lookup has duplicate IDs or P-codes")
                existing[pcode] = raster_id
                historical_rows[pcode] = row

    next_id = max(existing.values(), default=0) + 1
    current_by_pcode = {record["adm3_pcode"]: record for record in records}
    for pcode in sorted(current_by_pcode):
        if pcode not in existing:
            existing[pcode] = next_id
            next_id += 1
    if max(existing.values(), default=0) > 65535:
        raise RuntimeError("woreda raster IDs exceed uint16 logical capacity")

    rows = []
    for pcode, raster_id in existing.items():
        if pcode in current_by_pcode:
            record = current_by_pcode[pcode]
            row = {field: record.get(field, "") for field in LOOKUP_FIELDS}
            row.update(
                {"raster_id": raster_id, "adm3_pcode": pcode, "active": 1}
            )
            record["raster_id"] = raster_id
        else:
            row = dict(historical_rows[pcode])
            row["active"] = 0
        rows.append(row)
    rows.sort(key=lambda row: int(row["raster_id"]))
    return rows


def make_memory_woreda_layer(records: list[dict]):
    dataset = ogr.GetDriverByName("Memory").CreateDataSource("")
    layer = dataset.CreateLayer(
        "woredas", srs=expected_spatial_ref(), geom_type=ogr.wkbUnknown
    )
    layer.CreateField(ogr.FieldDefn("raster_id", ogr.OFTInteger))
    definition = layer.GetLayerDefn()
    for record in sorted(records, key=lambda item: item["adm3_pcode"]):
        feature = ogr.Feature(definition)
        feature.SetField("raster_id", record["raster_id"])
        feature.SetGeometry(record["geometry"])
        layer.CreateFeature(feature)
        feature = None
    return dataset, layer


def memory_raster(data_type, fill_value=0):
    dataset = gdal.GetDriverByName("MEM").Create("", COLS, ROWS, 1, data_type)
    dataset.SetGeoTransform(TRANSFORM.to_gdal())
    dataset.SetSpatialRef(expected_spatial_ref())
    dataset.SetMetadataItem("AREA_OR_POINT", "Area")
    dataset.GetRasterBand(1).Fill(fill_value)
    return dataset


def create_output(path: Path, *, description: str, nodata: float | None):
    remove_own_output(path)
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path),
        COLS,
        ROWS,
        1,
        gdal.GDT_Float32,
        options=(
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3",
        ),
    )
    dataset.SetGeoTransform(TRANSFORM.to_gdal())
    dataset.SetSpatialRef(expected_spatial_ref())
    dataset.SetMetadataItem("AREA_OR_POINT", "Area")
    band = dataset.GetRasterBand(1)
    band.SetDescription(description)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    return dataset


def write_lookup(rows: list[dict]) -> None:
    remove_own_output(LOOKUP_TEMP)
    with LOOKUP_TEMP.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOOKUP_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    before = {
        str(path.relative_to(ROOT)): source_snapshot(path)
        for path in (ADMIN0_PATH, WOREDA_PATH)
    }
    admin_dataset, admin_layer = open_layer(ADMIN0_PATH)
    woreda_dataset, woreda_layer = open_layer(WOREDA_PATH)
    if admin_layer.GetFeatureCount() != 1:
        raise RuntimeError("admin0 source must contain exactly one feature")

    records = read_woredas(woreda_layer)
    lookup_rows = assign_ids(records)
    memory_woredas, memory_woreda_layer = make_memory_woreda_layer(records)

    mask_memory = memory_raster(gdal.GDT_Byte, 0)
    admin_layer.ResetReading()
    gdal.RasterizeLayer(mask_memory, (1,), admin_layer, burn_values=(1,))
    mask = mask_memory.GetRasterBand(1).ReadAsArray().astype(bool)

    overlap_memory = memory_raster(gdal.GDT_UInt16, 0)
    memory_woreda_layer.ResetReading()
    gdal.RasterizeLayer(
        overlap_memory,
        (1,),
        memory_woreda_layer,
        burn_values=(1,),
        options=("MERGE_ALG=ADD",),
    )
    coverage_count = overlap_memory.GetRasterBand(1).ReadAsArray()
    overlap_cells = int((coverage_count > 1).sum())
    if overlap_cells:
        raise RuntimeError(
            f"woreda polygons overlap at {overlap_cells} canonical cell centres"
        )

    mask_output = create_output(
        MASK_TEMP, description="ethiopia_mask", nodata=None
    )
    mask_output.GetRasterBand(1).SetMetadataItem("logical_dtype", "uint8")
    mask_output.GetRasterBand(1).WriteArray(mask.astype("float32"))

    woreda_output = create_output(
        WOREDA_TEMP, description="woreda_id", nodata=NODATA
    )
    woreda_output.GetRasterBand(1).SetMetadataItem("logical_dtype", "uint16")
    woreda_output.GetRasterBand(1).Fill(NODATA)
    memory_woreda_layer.ResetReading()
    gdal.RasterizeLayer(
        woreda_output,
        (1,),
        memory_woreda_layer,
        options=("ATTRIBUTE=raster_id",),
    )
    woreda = woreda_output.GetRasterBand(1).ReadAsArray()
    uncovered_inside = mask & (woreda == NODATA)
    woreda[uncovered_inside] = 0
    woreda[~mask] = NODATA
    woreda_output.GetRasterBand(1).WriteArray(woreda)

    mask_output.FlushCache()
    woreda_output.FlushCache()
    mask_output = None
    woreda_output = None
    write_lookup(lookup_rows)

    assert_matches_grid(MASK_TEMP)
    assert_matches_grid(WOREDA_TEMP)
    remove_own_output(MASK_PATH)
    remove_own_output(WOREDA_ID_PATH)
    remove_own_output(LOOKUP_PATH)
    MASK_TEMP.replace(MASK_PATH)
    WOREDA_TEMP.replace(WOREDA_ID_PATH)
    LOOKUP_TEMP.replace(LOOKUP_PATH)

    legal_ids = {int(row["raster_id"]) for row in lookup_rows if row["active"] == 1}
    observed_ids = {int(value) for value in np.unique(woreda[mask]) if value > 0}
    if not observed_ids <= legal_ids:
        raise RuntimeError("woreda raster contains IDs absent from the lookup")
    zero_cell_woredas = [
        {
            "raster_id": int(row["raster_id"]),
            "adm3_pcode": row["adm3_pcode"],
            "adm3_name": row["adm3_name"],
        }
        for row in lookup_rows
        if int(row["active"]) and int(row["raster_id"]) not in observed_ids
    ]

    admin_dataset = woreda_dataset = memory_woredas = None
    mask_memory = overlap_memory = None
    after = {
        str(path.relative_to(ROOT)): source_snapshot(path)
        for path in (ADMIN0_PATH, WOREDA_PATH)
    }
    if before != after:
        raise RuntimeError("administrative source size or modification time changed")

    report = {
        "status": "passed",
        "sources_unchanged": True,
        "owner_approval": (
            "v04/2025-01-01 polygons approved for full-resolution "
            "lookup/rasterization and publication"
        ),
        "mask": {
            "path": str(MASK_PATH.relative_to(ROOT)),
            "inside_cells": int(mask.sum()),
            "outside_cells": int(mask.size - mask.sum()),
            "logical_dtype": "uint8",
            "physical_dtype": "float32",
            "nodata": None,
        },
        "woreda_id": {
            "path": str(WOREDA_ID_PATH.relative_to(ROOT)),
            "source_features": len(records),
            "active_lookup_rows": sum(
                int(row["active"]) for row in lookup_rows
            ),
            "historical_lookup_rows": sum(
                not int(row["active"]) for row in lookup_rows
            ),
            "observed_positive_ids": len(observed_ids),
            "zero_cell_woredas": zero_cell_woredas,
            "uncovered_inside_cells_set_to_zero": int(uncovered_inside.sum()),
            "overlap_cells": overlap_cells,
            "logical_dtype": "uint16",
            "physical_dtype": "float32",
            "outside_nodata": NODATA,
        },
        "lookup": str(LOOKUP_PATH.relative_to(ROOT)),
        "rasterization": "cell centre (GDAL default; ALL_TOUCHED disabled)",
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
