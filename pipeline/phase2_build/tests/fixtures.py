"""Small synthetic contracts shared by header and raster-build tests."""

import copy
from dataclasses import replace
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

from config import stack, stack_v2
from scoring.config import Crop
from scoring.engine import ScoringEngine
from scoring.soil import calibrate_selected_cv
from scoring.tests.fixtures import config_fixture
from scoring.tests.fixtures import values_fixture


def v1_manifest_fixture():
    return {"stack_version": stack.STACK_VERSION, "grid": stack.grid_metadata(),
            "overview_factors": stack.overview_factors(stack.grid.COLS, stack.grid.ROWS),
            "raster_file": "synthetic-v1.tif", "raster_sha256": "1" * 64,
            "bands": [{"index": index, "name": name, "description": f"Synthetic {name}",
                       "unit": stack.CORE_UNITS[name], "physical_dtype": "float32",
                       "logical_dtype": stack.logical_dtype(name), "nodata": stack.grid.NODATA,
                       "overview_resample": stack.overview_method(name),
                       "source": "AfricanSoilGrids" if name == "sg_p" else "Synthetic test source",
                       "source_file": "synthetic.tif", "resample": "none", "transform_applied": "none",
                       "missing_value_semantics": "Synthetic fixture only"}
                      for index, name in enumerate(stack.CORE_BAND_NAMES, start=1)]}


def engine_fixture(root: Path, *, six_crops=False):
    config = config_fixture(root)
    config = replace(config, band_units=dict(stack.CORE_UNITS))
    if six_crops:
        names = ("teff", "wheat", "maize", "sorghum", "barley", "chickpea")
        config = replace(config, crops=tuple(Crop(name, name.title(), ("kiremt", "belg"), index) for index, name in enumerate(names, start=1)),
                         requirements={name: copy.deepcopy(config.requirements["wheat"]) for name in names})
    calibration = calibrate_selected_cv(np.array([0.1, 0.2, 0.3, 0.8], dtype=np.float32), missing_cells=0,
                                        low_percentile=25, high_percentile=90, method="linear")
    return ScoringEngine(config, calibration)


def v2_manifest_fixture(root: Path, *, six_crops=True):
    engine = engine_fixture(root, six_crops=six_crops)
    v1 = v1_manifest_fixture()
    return {"stack_version": stack_v2.STACK_VERSION, "grid": stack.grid_metadata(),
            "overview_factors": v1["overview_factors"], "raster_file": "synthetic-v2.tif",
            "raster_sha256": "3" * 64, "raster_size_bytes": 1, "build_script": "pipeline/phase2_build/build_scores.py",
            "parent_stack": {"stack_version": 1, "raster_file": v1["raster_file"], "raster_sha256": v1["raster_sha256"],
                             "manifest_file": "synthetic-v1.json", "manifest_sha256": "2" * 64,
                             "band_records_sha256": stack_v2.metadata_sha256(v1["bands"])},
            "scoring": stack_v2.scoring_metadata(engine),
            "build_inputs": {row["file"]: {"sha256": row["sha256"], "size_bytes": Path(row["file"]).stat().st_size} for row in engine.config.inputs.values()},
            "bands": copy.deepcopy(v1["bands"]) + stack_v2.score_records(engine, v1["raster_file"])}


def source_fixture(width=11, height=9):
    """Invented tiny data, including gaps/exclusions; no production files read."""
    source = gdal.GetDriverByName("MEM").Create("", width, height, len(stack.CORE_BAND_NAMES), gdal.GDT_Float32)
    source.SetGeoTransform(stack.grid.TRANSFORM.to_gdal())
    crs = osr.SpatialReference()
    crs.SetFromUserInput(stack.grid.CRS)
    source.SetSpatialRef(crs)
    source.SetMetadataItem("AREA_OR_POINT", "Area")
    source.SetMetadataItem("stack_version", "1")
    bands = {name: np.full((height, width), index / 10, dtype=np.float32) for index, name in enumerate(stack.CORE_BAND_NAMES, start=1)}
    for name, value in values_fixture().items():
        bands[name].fill(value)
    bands["woreda_id"].fill(1)
    bands["ethiopia_mask"][0, :] = 0
    bands["ethiopia_mask"][:, 0] = 0
    bands["kiremt_cv"][:] = np.linspace(0.1, 0.7, width, dtype=np.float32)[None, :]
    bands["es_confidence"][:, 2::3] = 0
    bands["es_confidence"][2, 4] = 0.25
    for nutrient in stack.ETHIOSIS_NUTRIENTS:
        bands[f"es_{nutrient}"][bands["es_confidence"] == 0] = stack.grid.NODATA
    bands["landcover"][1, 1] = 50
    bands["landcover"][1, 2] = 80
    bands["sg_phh2o"][1, 1] = stack.grid.NODATA
    bands["sg_phh2o"][1, 3] = stack.grid.NODATA
    bands["landcover"][1, 4] = stack.grid.NODATA
    bands["slope_mean_deg"][1, 5] = 45
    bands["kiremt_mean_mm"][1, 7] = 150
    bands["belg_mean_mm"][1, 7] = 800
    bands["elevation_m"][2, 1] = 3200
    bands["elevation_m"][2, 2] = 800
    bands["elevation_m"][2, 3] = 100
    inside = bands["ethiopia_mask"] == 1
    records = v1_manifest_fixture()["bands"]
    for record in records:
        name = record["name"]
        if name != "ethiopia_mask":
            bands[name][~inside] = stack.grid.NODATA
        band = source.GetRasterBand(record["index"])
        band.SetDescription(name)
        band.SetUnitType(record["unit"])
        band.SetNoDataValue(stack.grid.NODATA)
        for key in ("unit", "source", "source_file", "physical_dtype", "logical_dtype", "overview_resample"):
            band.SetMetadataItem(key, record[key])
        band.SetMetadataItem("opaque_original_provenance", "Preserve this original v1 tag")
        band.WriteArray(bands[name])
    return source


def calibrated_fixture(root: Path, source, *, six_crops=True):
    from pipeline.phase2_build.raster import calibrate_stack

    manifest = v2_manifest_fixture(root, six_crops=six_crops)
    frozen = stack_v2.frozen_engine(manifest)
    engine = ScoringEngine(frozen.config, calibrate_stack(source, stack.band_indices_by_name(source), frozen.config.params, window_size=3))
    manifest["scoring"] = stack_v2.scoring_metadata(engine)
    manifest["bands"] = manifest["bands"][:len(stack.CORE_BAND_NAMES)] + stack_v2.score_records(engine, manifest["parent_stack"]["raster_file"])
    return engine, manifest
