"""Small synthetic contracts shared by header and raster-build tests."""

import copy
from dataclasses import replace
from pathlib import Path

import numpy as np

from config import stack, stack_v2
from scoring.config import Crop
from scoring.engine import ScoringEngine
from scoring.soil import calibrate_selected_cv
from scoring.tests.fixtures import config_fixture


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
