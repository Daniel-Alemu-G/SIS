"""Read-only acceptance of all pixels, overviews and scalar/window equality."""

from __future__ import annotations

import numpy as np
from osgeo_utils.samples.validate_cloud_optimized_geotiff import validate as validate_cog

from config import grid, stack, stack_v2
from pipeline.phase2_build.diagnostics import band_statistics, collect_band
from pipeline.phase2_build.raster import bitwise_mismatches, calibrate_stack, encode_score, read_named, validate_physical, verify_overviews, windows
from scoring.engine import decode

# Engineering-only reproducibility/sample settings, not scientific tuning.
VERIFY_WINDOW_SIZE = 127
RANDOM_INSIDE_CELLS = 32
RANDOM_SEED = 0


def compare_pixels(source, candidate, manifest: dict, *, raw=None, window_size=VERIFY_WINDOW_SIZE) -> tuple[list[dict], dict]:
    engine = stack_v2.frozen_engine(manifest)
    parent_indices, indices = stack.band_indices_by_name(source), stack.band_indices_by_name(candidate)
    raw_indices = stack.band_indices_by_name(raw) if raw is not None else {}
    expected_names = (*stack.CORE_BAND_NAMES, *stack_v2.score_band_names(engine.config.crops))
    if tuple(indices) != expected_names or tuple(parent_indices) != stack.CORE_BAND_NAMES:
        raise ValueError("pixel validation requires the exact v1/v2 named-band order")
    if (candidate.RasterXSize, candidate.RasterYSize) != (source.RasterXSize, source.RasterYSize):
        raise ValueError("candidate dimensions differ from the immutable v1 source")
    if calibrate_stack(source, parent_indices, engine.config.params, window_size=window_size) != engine.calibration:
        raise ValueError("frozen CV calibration differs from fresh input-window calibration")
    legal = {name: np.unique(collect_band(source, parent_indices, name)) for name in stack.MODE_BANDS}
    totals = {name: {"name": name, "v1_copy_mismatches": 0 if name in parent_indices else None,
                     "shared_engine_mismatches": 0 if name not in parent_indices else None,
                     "raw_copy_mismatches": 0 if raw is not None else None,
                     "valid_inside_cells": 0, "missing_inside_cells": 0} for name in indices}
    representatives = {}
    matched = {}
    inside_count = 0
    width = source.RasterXSize
    for window in windows(width, source.RasterYSize, window_size):
        col, row, ncols, nrows = window
        values = read_named(source, parent_indices, stack.CORE_BAND_NAMES, window)
        actual = read_named(candidate, indices, expected_names, window)
        computed = engine.evaluate(values)
        raw_values = read_named(raw, raw_indices, expected_names, window) if raw is not None else {}
        inside = values["ethiopia_mask"] == 1
        inside_count += int(np.count_nonzero(inside))
        for record in manifest["bands"]:
            name = record["name"]
            if name in values:
                validate_physical(values[name], record, inside=inside, legal_values=legal.get(name))
            validate_physical(actual[name], record, inside=inside, legal_values=legal.get(name))
            expected = values[name] if name in values else encode_score(name, computed.bands[name])
            mismatches = bitwise_mismatches(actual[name], expected)
            if mismatches:
                raise ValueError(f"{name}: {mismatches} bitwise v1/shared-engine pixel mismatches")
            if raw is not None and bitwise_mismatches(actual[name], raw_values[name]):
                raise ValueError(f"{name}: raw assembly and COG pixels differ")
            valid = actual[name] != grid.NODATA
            totals[name]["valid_inside_cells"] += int(np.count_nonzero(inside & valid))
            totals[name]["missing_inside_cells"] += int(np.count_nonzero(inside & ~valid))
        weight = computed.bands["data_weight"]
        cutoff = engine.config.params["blend"]["detailed_label_min_weight"]
        basic_missing = np.zeros(inside.shape, dtype=bool)
        for factor in engine.config.basic_factors:
            basic_missing |= ~np.isfinite(decode(values[factor.input_band]))
        available_es = sum(np.isfinite(decode(values[factor.input_band])).astype(np.int16) for factor in engine.config.es_factors)
        slope_rule = engine.generic["slope"]
        conditions = {"outside_Ethiopia": ~inside, "basic_only": inside & (weight == 0),
                      "partial_support": inside & (weight > 0) & (weight < cutoff),
                      "detailed_support": inside & (weight >= cutoff),
                      "missing_basic_input": inside & basic_missing,
                      "partial_EthioSIS_nutrients": inside & (values["es_confidence"] > 0) & (available_es > 0) & (available_es < len(engine.config.es_factors)),
                      "built_up": inside & (values["landcover"] == 50), "open_water": inside & (values["landcover"] == 80),
                      "unknown_landcover": inside & (values["landcover"] == grid.NODATA),
                      "hard_zero_slope": inside & (slope_rule.score(decode(values[slope_rule.input_band])) == 0),
                      "no_best_crop_above_minimum": inside & (computed.bands["best_crop"] == 0),
                      "exact_best_crop_tie": inside & (computed.bands["best_crop_margin"] == 0),
                      "main_season_kiremt": inside & (computed.main_season == 1),
                      "main_season_belg": inside & (computed.main_season == 2)}
        for name, condition in conditions.items():
            count = int(np.count_nonzero(condition))
            matched[name] = matched.get(name, 0) + count
            if count:
                local_row, local_col = np.unravel_index(np.flatnonzero(condition)[0], condition.shape)
                flat = int((row + local_row) * width + col + local_col)
                representatives[name] = min(representatives.get(name, flat), flat)
    if inside_count != engine.calibration["valid_calibration_cells"] + engine.calibration["missing_calibration_cells"]:
        raise ValueError("calibration cell counts differ from the country mask")
    for name in stack.CORE_BAND_NAMES:
        original, copied = source.GetRasterBand(parent_indices[name]), candidate.GetRasterBand(indices[name])
        if original.GetMetadata() != copied.GetMetadata():
            raise ValueError(f"{name}: original v1 per-band metadata changed")
    return list(totals.values()), {"matched_cells": matched, "representatives": representatives,
                                  "inside_cells": inside_count, "outside_cells": width * source.RasterYSize - inside_count}


def scalar_checks(source, candidate, manifest: dict, representatives: dict, *, random_count=RANDOM_INSIDE_CELLS, seed=RANDOM_SEED) -> dict:
    engine = stack_v2.frozen_engine(manifest)
    parent_indices, indices = stack.band_indices_by_name(source), stack.band_indices_by_name(candidate)
    inside = collect_band(source, parent_indices, "ethiopia_mask") == 1
    choices = np.flatnonzero(inside)
    selected = {}
    for name, flat in representatives.items():
        selected.setdefault(flat, []).append(name)
    for flat in np.random.default_rng(seed).choice(choices, size=min(random_count, choices.size), replace=False):
        selected.setdefault(int(flat), []).append("seeded_random_inside")
    details = []
    for flat, selection in sorted(selected.items()):
        row, col = divmod(flat, source.RasterXSize)
        window = (col, row, 1, 1)
        scalar_values = {name: value[0, 0] for name, value in read_named(source, parent_indices, stack.CORE_BAND_NAMES, window).items()}
        scalar = engine.score_cell(scalar_values)
        stored = read_named(candidate, indices, stack_v2.score_band_names(engine.config.crops), window)
        for name, value in scalar["scores"].items():
            expected = np.array([[grid.NODATA if value is None else value]], dtype=np.float32)
            if bitwise_mismatches(stored[name], expected):
                raise ValueError(f"cell ({row}, {col}) {name}: scalar/window/stack mismatch")
        details.append({"row": row, "col": col, "selection": selection, "scores": scalar["scores"],
                        "main_season": scalar["main_season"], "soil_tier": scalar["soil"]["tier"]})
    return {"compared_cells": len(details), "random_inside_cells": min(random_count, int(choices.size)), "seed": seed,
            "bitwise_mismatches": 0, "cells": details}


def validate_package(source, candidate, manifest: dict, v1_manifest: dict, *, raw=None) -> dict:
    stack.assert_stack_matches_manifest(source, v1_manifest)
    indices = stack_v2.assert_stack_matches_manifest(candidate, manifest, v1_manifest=v1_manifest)
    if raw is not None:
        stack_v2.assert_stack_matches_manifest(raw, manifest, v1_manifest=v1_manifest, require_cog=False)
    warnings, errors, _ = validate_cog(candidate, full_check=True)
    if errors:
        raise ValueError(f"full COG structural validation failed: {errors}")
    bands, samples = compare_pixels(source, candidate, manifest, raw=raw)
    overviews = verify_overviews(candidate, manifest["bands"], manifest["overview_factors"], raw=raw, parent=source)
    # compare_pixels has proved every assembly pixel equals the COG. Reading
    # its band-interleaved data avoids decoding all COG bands for each statistic.
    statistics_dataset = raw if raw is not None else candidate
    statistics_indices = stack.band_indices_by_name(statistics_dataset)
    statistics = {}
    inside = collect_band(statistics_dataset, statistics_indices, "ethiopia_mask") == 1
    for record in manifest["bands"][len(stack.CORE_BAND_NAMES):]:
        actual = band_statistics(collect_band(statistics_dataset, statistics_indices, record["name"]), inside, record)
        if record.get("statistics") is not None and record["statistics"] != actual:
            raise ValueError(f"{record['name']}: frozen manifest statistics differ from actual values")
        statistics[record["name"]] = actual
    return {"status": "passed", "band_count": len(indices), "v1_band_count": len(stack.CORE_BAND_NAMES),
            "v1_manifest_prefix": "identical", "grid": "passed", "physical_storage": "passed",
            "country_mask": "passed", "cog_full_structural_check": "passed", "cog_warnings": warnings,
            "inside_cells": samples["inside_cells"], "outside_cells": samples["outside_cells"],
            "all_full_resolution_pixels_per_band": candidate.RasterXSize * candidate.RasterYSize,
            "verification_window_size": VERIFY_WINDOW_SIZE, "bands": bands, "overviews": overviews,
            "representative_categories": samples["matched_cells"],
            "scalar_checks": scalar_checks(source, candidate, manifest, samples["representatives"]),
            "statistics_read_source": "bitwise_verified_assembly" if raw is not None else "candidate_cog",
            "appended_band_statistics": statistics}
