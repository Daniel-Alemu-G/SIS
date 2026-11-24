"""The approved Phase 1 stack/manifest contract, with name-based access.

This validates storage and provenance metadata, not scientific scores. Grid
values come only from config.grid. New bands or reordered names require a new,
explicit version contract rather than silently changing common_stack_v1.
"""

import math
import re
from os import PathLike

from config import grid
from config.grid import assert_matches_grid

STACK_VERSION = 1
BLOCK_SIZE = 256
ETHIOSIS_NUTRIENTS = ("ca", "cu", "fe", "mg", "mn", "mo", "k", "si", "na", "s", "zn")
CORE_BAND_NAMES = (
    "ethiopia_mask", "woreda_id",
    "sg_phh2o", "sg_sand", "sg_clay", "sg_silt", "sg_cec", "sg_soc",
    "sg_n", "sg_cfvo", "sg_p",
    *(f"es_{name}" for name in ETHIOSIS_NUTRIENTS), "es_confidence",
    "elevation_m", "slope_mean_deg", "landcover",
    "kiremt_mean_mm", "kiremt_std_mm", "kiremt_cv",
    "belg_mean_mm", "belg_std_mm", "belg_cv",
    "bega_mean_mm", "bega_std_mm", "bega_cv",
    "tmean_kiremt_c", "tmean_bega_c", "tmean_belg_c",
)
MODE_BANDS = frozenset({"ethiopia_mask", "woreda_id", "landcover"})
CORE_UNITS = {
    "ethiopia_mask": "0 or 1", "woreda_id": "id", "sg_phh2o": "pH",
    "sg_sand": "%", "sg_clay": "%", "sg_silt": "%",
    "sg_cec": "cmol(c)/kg", "sg_soc": "g/kg", "sg_n": "mg/kg",
    "sg_cfvo": "vol %", "sg_p": "mg/100kg",
    **{f"es_{name}": "mg/kg" for name in ETHIOSIS_NUTRIENTS},
    "es_confidence": "fraction", "elevation_m": "m",
    "slope_mean_deg": "degrees", "landcover": "class code",
    **{f"{season}_{stat}": unit for season in ("kiremt", "belg", "bega")
       for stat, unit in (("mean_mm", "mm"), ("std_mm", "mm"), ("cv", "fraction"))},
    **{f"tmean_{season}_c": "degrees Celsius" for season in ("kiremt", "bega", "belg")},
}


class StackContractError(ValueError):
    """Raised when a v1 raster or manifest violates its approved contract."""


def grid_metadata() -> dict:
    """Serialize the canonical grid without duplicating any numeric settings."""
    return {
        "crs": grid.CRS,
        "origin": [grid.WEST, grid.NORTH],
        "cell_deg": grid.RES,
        "cols": grid.COLS,
        "rows": grid.ROWS,
        "transform_gdal": list(grid.TRANSFORM.to_gdal()),
        "extent": {"west": grid.WEST, "north": grid.NORTH,
                   "east": grid.EAST, "south": grid.SOUTH},
        "pixel_convention": "pixel-is-area",
    }


def logical_dtype(name: str) -> str:
    if name in {"ethiopia_mask", "landcover"}:
        return "uint8"
    if name == "woreda_id":
        return "uint16"
    return "float32"


def overview_method(name: str) -> str:
    return "mode" if name in MODE_BANDS else "average"


def overview_factors(width: int, height: int) -> list[int]:
    """Use powers of two until the smallest overview fits one COG tile."""
    factors = []
    factor = 2
    while max(width, height) > BLOCK_SIZE:
        factors.append(factor)
        if math.ceil(max(width, height) / factor) <= BLOCK_SIZE:
            break
        factor *= 2
    return factors


def band_indices_by_name(dataset, *, unnamed_alias: str | None = None) -> dict[str, int]:
    """Resolve descriptions, allowing only an explicit unnamed-singleton alias."""
    if unnamed_alias is not None and dataset.RasterCount != 1:
        raise StackContractError("an unnamed source alias requires exactly one band")
    indices = {}
    for index in range(1, dataset.RasterCount + 1):
        name = dataset.GetRasterBand(index).GetDescription() or unnamed_alias
        if not name or name in indices:
            raise StackContractError(f"missing or duplicate raster band name: {name!r}")
        indices[name] = index
    return indices


def validate_manifest(manifest: dict) -> dict[str, dict]:
    """Validate the fixed v1 order, units, logical types, and physical encoding."""
    if not isinstance(manifest, dict):
        raise StackContractError("manifest must be an object")
    version = manifest.get("stack_version")
    if type(version) is not int or version != STACK_VERSION:
        raise StackContractError(f"manifest stack_version must be {STACK_VERSION}")
    if manifest.get("grid") != grid_metadata():
        raise StackContractError("manifest grid differs from pipeline/grid.yaml")
    expected_factors = overview_factors(grid.COLS, grid.ROWS)
    if manifest.get("overview_factors") != expected_factors:
        raise StackContractError("manifest overview factors differ from the v1 contract")
    records = manifest.get("bands")
    if not isinstance(records, list) or len(records) != len(CORE_BAND_NAMES):
        raise StackContractError(f"manifest must contain {len(CORE_BAND_NAMES)} core bands")
    result = {}
    for index, (record, expected_name) in enumerate(zip(records, CORE_BAND_NAMES), start=1):
        if not isinstance(record, dict):
            raise StackContractError(f"manifest band {index} must be an object")
        name = record.get("name")
        if name != expected_name or not re.fullmatch(r"[a-z][a-z0-9_]*", name):
            raise StackContractError(f"manifest band {index} must be {expected_name!r}")
        if type(record.get("index")) is not int or record["index"] != index:
            raise StackContractError(f"{name}: invalid manifest index")
        expected = {
            "unit": CORE_UNITS[name], "physical_dtype": "float32",
            "logical_dtype": logical_dtype(name), "nodata": grid.NODATA,
            "overview_resample": overview_method(name),
        }
        for key, value in expected.items():
            if record.get(key) != value:
                raise StackContractError(f"{name}: {key} must be {value!r}")
        for key in ("description", "source", "source_file", "resample", "missing_value_semantics"):
            value = record.get(key)
            if not isinstance(value, str) or not value.strip() or value.strip().upper() == "TBD":
                raise StackContractError(f"{name}: missing or unresolved {key}")
        if name == "sg_p":
            if "AfricanSoilGrids" not in record["source"]:
                raise StackContractError("sg_p must identify AfricanSoilGrids, not global SoilGrids")
            if record.get("transform_applied") != "none":
                raise StackContractError("sg_p must retain its owner-approved values without scaling")
        result[name] = record
    return result


def assert_stack_matches_manifest(dataset_or_path, manifest: dict, *, require_cog: bool = True) -> dict[str, int]:
    """Check the real grid, version, bands/tags, and COG header read-only.

    Pixel/source and overview-content comparisons are separate acceptance
    checks; a COG layout marker alone is not a complete raster validation.
    """
    from osgeo import gdal

    records = validate_manifest(manifest)
    opened_here = isinstance(dataset_or_path, (str, PathLike))
    if opened_here:
        dataset = gdal.OpenEx(str(dataset_or_path), gdal.OF_RASTER | gdal.OF_READONLY)
        if dataset is None:
            raise StackContractError(f"cannot open stack read-only: {dataset_or_path}")
    else:
        dataset = dataset_or_path
    assert_matches_grid(dataset)
    if dataset.GetMetadataItem("stack_version") != str(manifest["stack_version"]):
        raise StackContractError("raster stack_version does not match manifest")
    indices = band_indices_by_name(dataset)
    if tuple(indices) != CORE_BAND_NAMES:
        raise StackContractError("raster descriptions/order do not match the v1 manifest")
    if require_cog:
        image = dataset.GetMetadata("IMAGE_STRUCTURE")
        if image.get("LAYOUT") != "COG" or image.get("COMPRESSION") != "DEFLATE" or image.get("PREDICTOR") != "3":
            raise StackContractError("stack must be a DEFLATE/predictor-3 COG")
    for name, record in records.items():
        band = dataset.GetRasterBand(indices[name])
        if band.DataType != gdal.GDT_Float32 or band.GetNoDataValue() != grid.NODATA:
            raise StackContractError(f"{name}: physical dtype/nodata mismatch")
        if band.GetUnitType() != record["unit"]:
            raise StackContractError(f"{name}: raster unit does not match manifest")
        for key in ("unit", "source", "source_file", "logical_dtype", "overview_resample"):
            if band.GetMetadataItem(key) != record[key]:
                raise StackContractError(f"{name}: raster {key} does not match manifest")
        if require_cog:
            if band.GetBlockSize() != [BLOCK_SIZE, BLOCK_SIZE]:
                raise StackContractError(f"{name}: unexpected COG block size")
            if band.GetOverviewCount() != len(manifest["overview_factors"]):
                raise StackContractError(f"{name}: wrong overview count")
            for overview_index, factor in enumerate(manifest["overview_factors"]):
                overview = band.GetOverview(overview_index)
                if (overview.XSize, overview.YSize) != (math.ceil(grid.COLS / factor), math.ceil(grid.ROWS / factor)):
                    raise StackContractError(f"{name}: wrong overview dimensions")
    if opened_here:
        dataset = None
    return indices
