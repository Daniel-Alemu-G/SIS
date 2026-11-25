"""Explicit append-only scoring contract; the fixed v1 validator is unchanged.

Logical categorical nodata (255) never appears in the physical float32 stack:
it is encoded as the canonical -9999, like all other missing stack values.
The manifest freezes the scoring inputs/calibration needed by point reports.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
from os import PathLike

from config import grid, stack
from scoring.config import ConfigurationError, from_metadata
from scoring.engine import ScoringEngine

STACK_VERSION = 2
FORMULA_VERSION = "phase2-1"
SOIL_BAND_NAMES = ("soil_quality_basic", "soil_quality_detailed", "soil_quality", "soil_limiting_factor",
                   "data_weight", "drought_risk", "soil_score_overall")
BEST_CROP_BAND_NAMES = ("best_crop", "best_crop_score", "best_crop_margin")
TIES = {
    "soil_factor": "First lowest factor in active soil CSV order.",
    "crop_factor": "First lowest factor in gate order, then basic additional order, then active EthioSIS CSV order.",
    "main_season": "Higher kiremt/belg mean rainfall; an exact mean tie chooses kiremt; bega excluded.",
    "crop_season": "Highest rainfall score, then higher mean rainfall, then the crop CSV season order.",
    "crop": "Lowest original 1-based crops.csv index wins an exact viability tie; skipped crops are never renumbered.",
}
MISSING_DATA = {
    "basic_and_gates": "Any missing active basic/gate input propagates; even hard zeros do not conceal an unknown active input.",
    "ethiosis": "Ignore missing active nutrients in the detailed tier and scale confidence by their availability fraction; no available nutrient means missing detailed quality.",
    "landcover": "Open water (80), built-up (50) and unknown land cover have no crop viability; current cropland is not a filter.",
    "best_crop": "Ignore missing crops; all missing is logical 255/physical -9999; below minimum is valid 0; a sole valid crop has a missing runner-up margin.",
    "categorical": "Valid IDs are exact integers in float32. Logical 255 is encoded as physical -9999, never published as a class.",
    "outside": "Every band except the meaningful 0/1 country mask is physical -9999 outside Ethiopia; JSON missing values are null.",
}
RELIABILITY_MEANING = "Normal approximation to P(rain >= configured fraction of seasonal mean), in percent; NOT probability of meeting a crop water requirement. Raw seasonal CV is unchanged."


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def metadata_sha256(value) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def score_band_names(crops) -> tuple[str, ...]:
    return (*SOIL_BAND_NAMES, *(name for crop in crops for name in (f"viability_{crop.crop}", f"viability_{crop.crop}_limiting")),
            *BEST_CROP_BAND_NAMES)


def scoring_metadata(engine: ScoringEngine) -> dict:
    return {"formula_version": FORMULA_VERSION, "configuration": engine.config.metadata(),
            "calibration": copy.deepcopy(engine.calibration), "ties": dict(TIES), "missing_data": dict(MISSING_DATA),
            "reliability_meaning": RELIABILITY_MEANING}


def score_records(engine: ScoringEngine, source_file: str) -> list[dict]:
    """Append the seven soil bands, crop pairs in CSV order, and best-crop trio."""
    descriptions = {
        "soil_quality_basic": "Crop-independent soil quality from active SoilGrids and mean-slope factors",
        "soil_quality_detailed": "Basic soil factors plus available active EthioSIS IDW nutrient estimates",
        "soil_quality": "Smooth confidence/completeness blend of basic and detailed soil quality",
        "soil_limiting_factor": "Lowest factor in the dominant soil tier, using the frozen factor lookup",
        "data_weight": "EthioSIS proximity confidence times the fraction of available active soil nutrients",
        "drought_risk": "Variability of the larger-mean rainy season, calibrated by Ethiopian CV percentiles",
        "soil_score_overall": "Soil quality with rainfall-variability penalty; soil and rainfall reliability, not pure soil",
        "best_crop": "Original 1-based crop ID with highest viability; 0 means below the approved minimum",
        "best_crop_score": "Highest valid crop viability, including scores below the best-crop minimum",
        "best_crop_margin": "Highest minus second-highest valid crop viability; missing without a runner-up",
    }
    lookups = {"soil_limiting_factor": {str(index): name for name, index in engine.soil_lookup.items()},
               "best_crop": {"0": "none", **{str(crop.index): crop.crop for crop in engine.config.crops}}}
    for crop in engine.config.crops:
        descriptions[f"viability_{crop.crop}"] = f"{crop.display_name} advisory viability: hybrid hard-zero gate times fertility/terrain modifier"
        descriptions[f"viability_{crop.crop}_limiting"] = f"{crop.display_name} limiting gate/additional factor; hard-zero slope takes precedence"
        lookups[f"viability_{crop.crop}_limiting"] = {str(index): name for name, index in engine.crop_lookups[crop.crop].items()}
    records = []
    for index, name in enumerate(score_band_names(engine.config.crops), start=len(stack.CORE_BAND_NAMES) + 1):
        categorical = name in lookups
        missing = MISSING_DATA["best_crop"] if name in BEST_CROP_BAND_NAMES else MISSING_DATA["categorical"] if categorical else "Missing active basic/gate input or outside Ethiopia is physical -9999/computational NaN/API null; no imputation."
        if name == "soil_quality_detailed":
            missing += " Missing where no active EthioSIS nutrient is available."
        if name == "data_weight":
            missing = "0 is valid basic-only support inside Ethiopia; missing confidence is treated as zero. Outside Ethiopia is -9999."
        if name == "drought_risk" or name == "soil_score_overall":
            missing += " Missing main-season comparison means or selected CV also propagates."
        if name.startswith("viability_"):
            missing += " " + MISSING_DATA["landcover"]
        legal = sorted(int(key) for key in lookups[name]) if categorical else None
        record = {"index": index, "name": name, "description": descriptions[name], "group": "scoring",
                  "unit": "id" if categorical else "fraction", "physical_dtype": "float32",
                  "logical_dtype": "uint8" if categorical else "float32", "logical_nodata": 255 if categorical else None,
                  "nodata": grid.NODATA, "overview_resample": "mode" if categorical else "average",
                  "source": "Owner-approved scoring configuration and accepted stack v1",
                  "source_file": source_file, "resample": "none; score aligned cells offline",
                  "transform_applied": "Shared scoring engine; formulas/thresholds frozen in scoring metadata",
                  "missing_value_semantics": missing, "valid_min": min(legal) if categorical else 0,
                  "valid_max": max(legal) if categorical else 1, "valid_range_role": "Strict score/ID contract, not a raw-input clipping rule.",
                  "resolution_label": "Nominal 1 km canonical cell; advisory score from estimates, not an exact 90 m value or measured sample.",
                  "scoring_metadata": "scoring.configuration and scoring.calibration",
                  "formula_version": FORMULA_VERSION}
        if categorical:
            record.update({"lookup": lookups[name], "valid_values": legal})
        records.append(record)
    return records


def frozen_engine(manifest: dict) -> ScoringEngine:
    metadata = manifest.get("scoring")
    expected = {"formula_version", "configuration", "calibration", "ties", "missing_data", "reliability_meaning"}
    if not isinstance(metadata, dict) or set(metadata) != expected or metadata["formula_version"] != FORMULA_VERSION:
        raise stack.StackContractError("missing or incompatible frozen scoring metadata")
    if metadata["ties"] != TIES or metadata["missing_data"] != MISSING_DATA:
        raise stack.StackContractError("scoring tie/missing-data policies differ from the v2 contract")
    if metadata["reliability_meaning"] != RELIABILITY_MEANING:
        raise stack.StackContractError("rainfall reliability interpretation is missing")
    try:
        config = from_metadata(metadata["configuration"], manifest)
    except (ConfigurationError, KeyError, TypeError, ValueError) as exc:
        raise stack.StackContractError(f"invalid frozen scoring configuration: {exc}") from exc
    calibration = metadata["calibration"]
    required = {"cv_low", "cv_high", "low_percentile", "high_percentile", "percentile_method", "valid_calibration_cells", "missing_calibration_cells", "selection"}
    if not isinstance(calibration, dict) or set(calibration) != required:
        raise stack.StackContractError("missing or incomplete CV calibration")
    for key in ("cv_low", "cv_high", "low_percentile", "high_percentile"):
        if not isinstance(calibration[key], (int, float)) or isinstance(calibration[key], bool) or not math.isfinite(calibration[key]):
            raise stack.StackContractError("CV calibration values must be finite numbers")
    if not 0 <= calibration["cv_low"] < calibration["cv_high"]:
        raise stack.StackContractError("CV calibration cut-offs must be ordered and nonnegative")
    drought = config.params["drought"]
    if (calibration["low_percentile"], calibration["high_percentile"], calibration["percentile_method"]) != (drought["cv_low_percentile"], drought["cv_high_percentile"], drought["percentile_method"]):
        raise stack.StackContractError("CV calibration does not match approved percentile parameters")
    if type(calibration["valid_calibration_cells"]) is not int or calibration["valid_calibration_cells"] <= 0 or type(calibration["missing_calibration_cells"]) is not int or calibration["missing_calibration_cells"] < 0:
        raise stack.StackContractError("CV calibration cell counts are invalid")
    if calibration["selection"] != "CV of the larger kiremt/belg mean; exact mean ties choose kiremt; bega excluded":
        raise stack.StackContractError("CV calibration season-selection policy differs")
    engine = ScoringEngine(config, calibration)
    if canonical_json(metadata) != canonical_json(scoring_metadata(engine)):
        raise stack.StackContractError("frozen scoring metadata is inconsistent")
    return engine


def _fingerprint(value, label: str) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise stack.StackContractError(f"{label}: invalid SHA-256 fingerprint")


def validate_manifest(manifest: dict, *, v1_manifest: dict | None = None) -> dict[str, dict]:
    if not isinstance(manifest, dict) or type(manifest.get("stack_version")) is not int or manifest["stack_version"] != STACK_VERSION:
        raise stack.StackContractError("manifest stack_version must be 2")
    records = manifest.get("bands")
    if not isinstance(records, list) or len(records) < len(stack.CORE_BAND_NAMES):
        raise stack.StackContractError("v2 manifest is missing its v1 prefix")
    prefix = records[:len(stack.CORE_BAND_NAMES)]
    stack.validate_manifest({"stack_version": stack.STACK_VERSION, "grid": manifest.get("grid"),
                             "overview_factors": manifest.get("overview_factors"), "bands": prefix})
    parent = manifest.get("parent_stack")
    if not isinstance(parent, dict) or set(parent) != {"stack_version", "raster_file", "raster_sha256", "manifest_file", "manifest_sha256", "band_records_sha256"} or type(parent["stack_version"]) is not int or parent["stack_version"] != stack.STACK_VERSION:
        raise stack.StackContractError("v2 requires an immutable stack-v1 parent fingerprint")
    for key in ("raster_sha256", "manifest_sha256", "band_records_sha256"):
        _fingerprint(parent[key], f"parent {key}")
    if parent["band_records_sha256"] != metadata_sha256(prefix):
        raise stack.StackContractError("v1 band records changed within the v2 manifest")
    for key in ("raster_file", "manifest_file"):
        if not isinstance(parent[key], str) or not parent[key]:
            raise stack.StackContractError("parent file paths are required")
    if v1_manifest is not None:
        stack.validate_manifest(v1_manifest)
        if prefix != v1_manifest["bands"]:
            raise stack.StackContractError("first 38 band records must be identical to the original v1 manifest")
        if "raster_sha256" in v1_manifest and parent["raster_sha256"] != v1_manifest["raster_sha256"]:
            raise stack.StackContractError("parent raster fingerprint differs from v1")
    engine = frozen_engine(manifest)
    expected = score_records(engine, parent["raster_file"])
    if len(records) != len(prefix) + len(expected):
        raise stack.StackContractError("v2 band count differs from the frozen crop/score contract")
    for record, template in zip(records[len(prefix):], expected, strict=True):
        if not isinstance(record, dict) or set(record) - (set(template) | {"statistics"}):
            raise stack.StackContractError("unknown v2 score-band metadata")
        for key, value in template.items():
            if record.get(key) != value or (key in {"index", "logical_nodata"} and value is not None and type(record.get(key)) is not int):
                raise stack.StackContractError(f"{template['name']}: {key} differs from the v2 scoring contract")
    for key in ("raster_file", "build_script"):
        if not isinstance(manifest.get(key), str) or not manifest[key]:
            raise stack.StackContractError(f"missing {key}")
    _fingerprint(manifest.get("raster_sha256"), "raster")
    if type(manifest.get("raster_size_bytes")) is not int or manifest["raster_size_bytes"] <= 0:
        raise stack.StackContractError("raster_size_bytes must be a positive integer")
    inputs = manifest.get("build_inputs")
    if not isinstance(inputs, dict) or not inputs:
        raise stack.StackContractError("build-input fingerprints are required")
    for path, row in inputs.items():
        if not isinstance(path, str) or not path or not isinstance(row, dict) or set(row) != {"sha256", "size_bytes"} or type(row["size_bytes"]) is not int or row["size_bytes"] < 0:
            raise stack.StackContractError("invalid build-input fingerprint record")
        _fingerprint(row["sha256"], path)
    for row in engine.config.inputs.values():
        if inputs.get(row["file"], {}).get("sha256") != row["sha256"]:
            raise stack.StackContractError("frozen scoring input hash differs from build provenance")
    return {record["name"]: record for record in records}


def assert_stack_matches_manifest(dataset_or_path, manifest: dict, *, v1_manifest: dict | None = None, require_cog: bool = True) -> dict[str, int]:
    """Read-only version/grid/storage/tag validation; pixels are checked offline."""
    from osgeo import gdal

    records = validate_manifest(manifest, v1_manifest=v1_manifest)
    opened_here = isinstance(dataset_or_path, (str, PathLike))
    dataset = gdal.OpenEx(str(dataset_or_path), gdal.OF_RASTER | gdal.OF_READONLY) if opened_here else dataset_or_path
    if dataset is None:
        raise stack.StackContractError(f"cannot open stack read-only: {dataset_or_path}")
    grid.assert_matches_grid(dataset)
    if dataset.GetMetadataItem("stack_version") != str(STACK_VERSION):
        raise stack.StackContractError("raster stack_version does not match manifest")
    for key, expected in {"parent_stack_sha256": manifest["parent_stack"]["raster_sha256"],
                          "scoring_configuration_sha256": metadata_sha256(manifest["scoring"]["configuration"]),
                          "formula_version": FORMULA_VERSION,
                          "cv_low": str(manifest["scoring"]["calibration"]["cv_low"]),
                          "cv_high": str(manifest["scoring"]["calibration"]["cv_high"])}.items():
        if dataset.GetMetadataItem(key) != expected:
            raise stack.StackContractError(f"raster {key} does not match frozen manifest")
    indices = stack.band_indices_by_name(dataset)
    if tuple(indices) != tuple(records):
        raise stack.StackContractError("raster descriptions/order do not match the v2 manifest")
    if require_cog:
        image = dataset.GetMetadata("IMAGE_STRUCTURE")
        if image.get("LAYOUT") != "COG" or image.get("COMPRESSION") != "DEFLATE" or image.get("PREDICTOR") != "3":
            raise stack.StackContractError("stack must be a DEFLATE/predictor-3 COG")
    for name, record in records.items():
        band = dataset.GetRasterBand(indices[name])
        if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != grid.NODATA:
            raise stack.StackContractError(f"{name}: physical dtype/nodata mismatch")
        if band.GetUnitType() != record["unit"]:
            raise stack.StackContractError(f"{name}: raster unit does not match manifest")
        for key in ("unit", "source", "source_file", "logical_dtype", "overview_resample"):
            if band.GetMetadataItem(key) != record[key]:
                raise stack.StackContractError(f"{name}: raster {key} does not match manifest")
        if "lookup" in record and band.GetMetadataItem("lookup") != canonical_json(record["lookup"]):
            raise stack.StackContractError(f"{name}: raster factor/crop lookup differs from manifest")
        if band.GetScale() not in (None, 1.0) or band.GetOffset() not in (None, 0.0):
            raise stack.StackContractError(f"{name}: an unexpected scale/offset would alter stored values")
        if require_cog:
            if band.GetBlockSize() != [stack.BLOCK_SIZE, stack.BLOCK_SIZE] or band.GetOverviewCount() != len(manifest["overview_factors"]):
                raise stack.StackContractError(f"{name}: incorrect COG tile/overview layout")
            for overview_index, factor in enumerate(manifest["overview_factors"]):
                overview = band.GetOverview(overview_index)
                if (overview.XSize, overview.YSize) != (math.ceil(grid.COLS / factor), math.ceil(grid.ROWS / factor)):
                    raise stack.StackContractError(f"{name}: wrong overview dimensions")
    return indices
