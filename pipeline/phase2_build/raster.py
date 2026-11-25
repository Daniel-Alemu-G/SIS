"""Bounded, name-based raster I/O and deterministic COG assembly.

Only stack v1 is read. No native raster, interpolation, reprojection or seasonal
aggregation is performed. All grid values come from config.grid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from osgeo import gdal, osr

from config import grid, stack, stack_v2
from scoring.engine import ScoringEngine, decode
from scoring.soil import calibrate_selected_cv, main_rainy_season

gdal.UseExceptions()


def windows(width: int, height: int, size: int = stack.BLOCK_SIZE):
    if any(type(value) is not int or value <= 0 for value in (width, height, size)):
        raise ValueError("window dimensions and size must be positive integers")
    for row in range(0, height, size):
        for col in range(0, width, size):
            yield col, row, min(size, width - col), min(size, height - row)


def read_named(dataset, indices: dict[str, int], names, window=None) -> dict[str, np.ndarray]:
    names = tuple(names)
    if not names or len(names) != len(set(names)) or set(names) - set(indices):
        raise ValueError("requested raster bands must be present, nonempty and unique by name")
    if window is None:
        window = (0, 0, dataset.RasterXSize, dataset.RasterYSize)
    array = dataset.ReadAsArray(*window, band_list=[indices[name] for name in names])
    if array is None or array.dtype != np.float32:
        raise ValueError("stack reads must yield physical float32 values")
    if array.ndim == 2:
        array = array[None, ...]
    if array.ndim != 3 or array.shape[0] != len(names):
        raise ValueError("unexpected named-band read shape")
    return dict(zip(names, array, strict=True))


def encode_score(name: str, values) -> np.ndarray:
    array = np.asarray(values, dtype=np.float32)
    if np.any(np.isinf(array)):
        raise ValueError(f"{name}: infinite computed values cannot be published as nodata")
    missing = np.isnan(array) | ((array == 255) if name == "best_crop" else False)
    return np.where(missing, grid.NODATA, array).astype(np.float32)


def bitwise_mismatches(actual: np.ndarray, expected: np.ndarray) -> int:
    if actual.dtype != np.float32 or expected.dtype != np.float32 or actual.shape != expected.shape:
        raise ValueError("bitwise comparisons require identical float32 shapes")
    return int(np.count_nonzero(actual.view(np.uint32) != expected.view(np.uint32)))


def validate_physical(values: np.ndarray, record: dict, *, inside=None, legal_values=None) -> None:
    name = record["name"]
    if values.dtype != np.float32 or not np.all(np.isfinite(values)):
        raise ValueError(f"{name}: physical storage must be finite float32 with -9999 nodata")
    if name == "ethiopia_mask":
        if not np.all(np.isin(values, (0, 1))):
            raise ValueError("ethiopia_mask must contain only meaningful 0/1 pixels")
    elif inside is not None and np.any(values[~inside] != grid.NODATA):
        raise ValueError(f"{name}: non-nodata values outside Ethiopia")
    valid = values[values != grid.NODATA]
    if record["logical_dtype"] != "float32":
        if not np.all(valid == np.rint(valid)):
            raise ValueError(f"{name}: fractional categorical/ID values")
        allowed = record.get("valid_values", legal_values)
        if allowed is not None and not np.all(np.isin(valid, tuple(allowed))):
            raise ValueError(f"{name}: illegal categorical/ID values")
    elif record.get("formula_version") == stack_v2.FORMULA_VERSION:
        if np.any((valid < 0) | (valid > 1)):
            raise ValueError(f"{name}: score outside [0, 1]")


def calibrate_stack(dataset, indices: dict[str, int], params: dict, *, window_size=stack.BLOCK_SIZE) -> dict:
    names = ("ethiopia_mask", "kiremt_mean_mm", "belg_mean_mm", "kiremt_cv", "belg_cv")
    chunks = []
    missing = 0
    for window in windows(dataset.RasterXSize, dataset.RasterYSize, window_size):
        values = read_named(dataset, indices, names, window)
        if not np.all(np.isin(values["ethiopia_mask"], (0, 1))):
            raise ValueError("calibration country mask is not binary")
        inside = values["ethiopia_mask"] == 1
        _, cv = main_rainy_season(*(decode(values[name]) for name in names[1:]))
        available = inside & np.isfinite(cv)
        missing += int(np.count_nonzero(inside & ~available))
        if np.any(available):
            chunks.append(cv[available])
    sample = np.concatenate(chunks) if chunks else np.array([], dtype=np.float32)
    drought = params["drought"]
    return calibrate_selected_cv(sample, missing_cells=missing, low_percentile=drought["cv_low_percentile"],
                                 high_percentile=drought["cv_high_percentile"], method=drought["percentile_method"])


def configure_output(dataset, engine: ScoringEngine, parent_sha256: str) -> None:
    dataset.SetGeoTransform(grid.TRANSFORM.to_gdal())
    spatial_ref = osr.SpatialReference()
    spatial_ref.SetFromUserInput(grid.CRS)
    dataset.SetSpatialRef(spatial_ref)
    for key, value in {"AREA_OR_POINT": "Area", "stack_version": str(stack_v2.STACK_VERSION),
                       "role": "Offline precomputed advisory soil/crop scoring stack",
                       "parent_stack_sha256": parent_sha256, "formula_version": stack_v2.FORMULA_VERSION,
                       "scoring_configuration_sha256": stack_v2.metadata_sha256(engine.config.metadata()),
                       "cv_low": str(engine.calibration["cv_low"]), "cv_high": str(engine.calibration["cv_high"])}.items():
        dataset.SetMetadataItem(key, value)


def apply_overview_policy(dataset, records: list[dict], factors: list[int]) -> None:
    if not factors:
        return
    if dataset.BuildOverviews("AVERAGE", factors) != gdal.CE_None:
        raise ValueError("could not create common overview levels")
    indices = stack.band_indices_by_name(dataset)
    # Never cascade averages: nodata and nondivisible edge dimensions require
    # regenerating every level from the original full-resolution contributors.
    for record in records:
        band = dataset.GetRasterBand(indices[record["name"]])
        for index in range(band.GetOverviewCount()):
            if gdal.RegenerateOverview(band, band.GetOverview(index), record["overview_resample"].upper()) != gdal.CE_None:
                raise ValueError(f"{record['name']}: overview regeneration failed")


def write_stack(source, records: list[dict], engine: ScoringEngine, raw_path: Path, cog_path: Path, *, parent_sha256: str,
                overview_factors: list[int], window_size: int = stack.BLOCK_SIZE) -> None:
    """Caller validates the canonical v1 input and ownership of output paths."""
    expected_names = (*stack.CORE_BAND_NAMES, *stack_v2.score_band_names(engine.config.crops))
    if tuple(record["name"] for record in records) != expected_names:
        raise ValueError("assembly records do not match the append-only named-band contract")
    indices = stack.band_indices_by_name(source)
    if tuple(indices) != stack.CORE_BAND_NAMES:
        raise ValueError("assembly source must contain only the original v1 bands")
    raw_path, cog_path = Path(raw_path), Path(cog_path)
    if raw_path.resolve() == cog_path.resolve() or any(path.resolve() == Path(source.GetDescription()).resolve() for path in (raw_path, cog_path) if source.GetDescription()):
        raise ValueError("assembly cannot overwrite its source or write in place")
    output = gdal.GetDriverByName("GTiff").Create(str(raw_path), source.RasterXSize, source.RasterYSize, len(records), gdal.GDT_Float32,
        options=["TILED=YES", f"BLOCKXSIZE={stack.BLOCK_SIZE}", f"BLOCKYSIZE={stack.BLOCK_SIZE}",
                 "COMPRESS=DEFLATE", "PREDICTOR=3", "INTERLEAVE=BAND", "NUM_THREADS=1"])
    if output is None:
        raise ValueError("could not create candidate assembly raster")
    configure_output(output, engine, parent_sha256)
    output_indices = {}
    for record in records:
        name = record["name"]
        output_indices[name] = record["index"]
        band = output.GetRasterBand(output_indices[name])
        band.SetDescription(name)
        band.SetNoDataValue(grid.NODATA)
        band.SetUnitType(record["unit"])
        if name in indices:
            original = source.GetRasterBand(indices[name])
            band.SetMetadata(original.GetMetadata())
        else:
            for key in ("unit", "source", "source_file", "logical_dtype", "physical_dtype", "overview_resample", "resample",
                        "transform_applied", "missing_value_semantics", "resolution_label", "formula_version"):
                band.SetMetadataItem(key, str(record[key]))
            if "lookup" in record:
                band.SetMetadataItem("lookup", stack_v2.canonical_json(record["lookup"]))
                band.SetMetadataItem("logical_nodata", str(record["logical_nodata"]))
    for window in windows(source.RasterXSize, source.RasterYSize, window_size):
        values = read_named(source, indices, stack.CORE_BAND_NAMES, window)
        computed = engine.evaluate(values).bands
        if set(computed) != set(stack_v2.score_band_names(engine.config.crops)):
            raise ValueError("scoring engine output names differ from the appended-band contract")
        inside = values["ethiopia_mask"] == 1
        for record in records:
            name = record["name"]
            pixels = values[name] if name in values else encode_score(name, computed[name])
            validate_physical(pixels, record, inside=inside)
            output.GetRasterBand(output_indices[name]).WriteArray(pixels, window[0], window[1])
    output.FlushCache()
    apply_overview_policy(output, records, overview_factors)
    output.FlushCache()
    output = None
    temporary = cog_path.with_name(cog_path.stem + ".tmp.tif")
    translated = gdal.Translate(str(temporary), str(raw_path), format="COG", creationOptions=[
        "COMPRESS=DEFLATE", "PREDICTOR=FLOATING_POINT", "OVERVIEW_COMPRESS=DEFLATE",
        "OVERVIEW_PREDICTOR=FLOATING_POINT", f"BLOCKSIZE={stack.BLOCK_SIZE}",
        "OVERVIEWS=FORCE_USE_EXISTING", "NUM_THREADS=1"])
    if translated is None:
        raise ValueError("COG translation failed")
    translated = None
    temporary.replace(cog_path)


def verify_overviews(dataset, records: list[dict], factors: list[int], *, raw=None, parent=None) -> list[dict]:
    """Compare every overview pixel to an independent, direct full-res reference."""
    indices = stack.band_indices_by_name(dataset)
    raw_indices = stack.band_indices_by_name(raw) if raw is not None else {}
    parent_indices = stack.band_indices_by_name(parent) if parent is not None else {}
    results = []
    for record in records:
        name = record["name"]
        band = dataset.GetRasterBand(indices[name])
        full = band.ReadAsArray()
        legal = np.unique(full[full != grid.NODATA]) if record["logical_dtype"] != "float32" else None
        reference = gdal.GetDriverByName("MEM").Create("", dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Float32)
        reference_band = reference.GetRasterBand(1)
        reference_band.SetNoDataValue(grid.NODATA)
        reference_band.WriteArray(full)
        levels = []
        for index, factor in enumerate(factors):
            actual_band = band.GetOverview(index)
            if actual_band is None:
                raise ValueError(f"{name}: missing overview {factor}")
            target = gdal.GetDriverByName("MEM").Create("", actual_band.XSize, actual_band.YSize, 1, gdal.GDT_Float32)
            target_band = target.GetRasterBand(1)
            target_band.SetNoDataValue(grid.NODATA)
            if gdal.RegenerateOverview(reference_band, target_band, record["overview_resample"].upper()) != gdal.CE_None:
                raise ValueError(f"{name}: independent overview generation failed")
            actual = actual_band.ReadAsArray()
            mismatch = bitwise_mismatches(actual, target_band.ReadAsArray())
            copied = bitwise_mismatches(actual, raw.GetRasterBand(raw_indices[name]).GetOverview(index).ReadAsArray()) if raw is not None else None
            inherited = bitwise_mismatches(actual, parent.GetRasterBand(parent_indices[name]).GetOverview(index).ReadAsArray()) if name in parent_indices else None
            if mismatch or copied or inherited:
                raise ValueError(f"{name} overview {factor}: reference={mismatch}, COG copy={copied}, v1={inherited} mismatches")
            validate_physical(actual, record, legal_values=legal)
            levels.append({"factor": factor, "shape_rows_cols": list(actual.shape), "method": record["overview_resample"],
                           "reference_mismatches": mismatch, "cog_copy_mismatches": copied, "v1_overview_mismatches": inherited})
            target = None
        reference = None
        results.append({"name": name, "levels": levels})
    return results
