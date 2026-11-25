"""Offline review evidence; diagnostic cut-offs are never scoring defaults.

Geographic elevation tests (3,000/1,200 m) and correlation flag (>0.95) are
specified review questions, not scientific thresholds introduced by this code.
Sensitivity values and score-class cut-offs come from the approved YAML.
"""

from __future__ import annotations

import copy
import math
from dataclasses import replace
from itertools import combinations

import numpy as np
from osgeo import gdal
from scipy.stats import rankdata, spearmanr

from config import grid, stack, stack_v2
from pipeline.phase2_build.raster import read_named, windows
from scoring import crops
from scoring.config import validate_params
from scoring.engine import ScoringEngine, decode

# Floating-point comparison tolerance only, not an agronomic threshold.
ROUNDING_TOLERANCE = 4 * np.finfo(np.float32).eps


def collect_band(dataset, indices: dict[str, int], name: str) -> np.ndarray:
    result = np.empty((dataset.RasterYSize, dataset.RasterXSize), dtype=np.float32)
    for window in windows(dataset.RasterXSize, dataset.RasterYSize):
        col, row, width, height = window
        result[row:row + height, col:col + width] = read_named(dataset, indices, (name,), window)[name]
    return result


def numeric_summary(values) -> dict:
    sample = np.asarray(values)
    sample = sample[np.isfinite(sample)]
    return {"count": int(sample.size), "min": float(sample.min()) if sample.size else None,
            "median": float(np.median(sample)) if sample.size else None,
            "mean": float(sample.mean(dtype=np.float64)) if sample.size else None,
            "p90": float(np.percentile(sample, 90, method="linear")) if sample.size else None,
            "max": float(sample.max()) if sample.size else None}


def band_statistics(values: np.ndarray, inside: np.ndarray, record: dict) -> dict:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{record['name']}: nonfinite physical values in statistics")
    valid = values != grid.NODATA
    domain = np.ones(inside.shape, dtype=bool) if record["name"] == "ethiopia_mask" else inside
    sample = values[domain & valid]
    stats = numeric_summary(sample)
    if record["logical_dtype"] != "float32":
        unique, counts = np.unique(sample, return_counts=True)
        histogram = {"kind": "categorical", "values": [int(value) for value in unique], "counts": [int(value) for value in counts]}
    elif not sample.size:
        histogram = {"kind": "empty", "edges": [], "counts": []}
    elif sample.min() == sample.max():
        histogram = {"kind": "constant", "values": [float(sample[0])], "counts": [int(sample.size)]}
    else:
        low, high = (0, 1) if record.get("formula_version") == stack_v2.FORMULA_VERSION else (float(sample.min()), float(sample.max()))
        counts, edges = np.histogram(sample, bins=np.linspace(low, high, 11))
        histogram = {"kind": "ten equal-width diagnostic bins", "edges": [float(value) for value in edges], "counts": [int(value) for value in counts]}
    return {**stats, "statistics_domain": "whole grid" if record["name"] == "ethiopia_mask" else "inside Ethiopia",
            "valid_inside_cells": int(np.count_nonzero(inside & valid)), "missing_inside_cells": int(np.count_nonzero(inside & ~valid)),
            "exact_zero_cells": int(np.count_nonzero(sample == 0)), "exact_one_cells": int(np.count_nonzero(sample == 1)),
            "below_sanity_range_cells": None if record.get("valid_min") is None else int(np.count_nonzero(sample < record["valid_min"])),
            "above_sanity_range_cells": None if record.get("valid_max") is None else int(np.count_nonzero(sample > record["valid_max"])),
            "histogram": histogram}


def score_classes(values, cutoffs_percent) -> np.ndarray:
    """An exact float32 cutoff enters the higher class; missing is logical 255."""
    values = decode(values)
    levels = np.asarray(cutoffs_percent, dtype=np.float64)
    if levels.ndim != 1 or not levels.size or not np.all(np.isfinite(levels)) or levels[0] <= 0 or levels[-1] >= 100 or np.any(np.diff(levels) <= 0):
        raise ValueError("class cutoffs must be finite, increasing percentages between 0 and 100")
    if np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
        raise ValueError("class inputs must be scores in [0, 1] or missing")
    cutoffs = (levels / 100).astype(np.float32)
    classes = np.searchsorted(cutoffs, values, side="right").astype(np.uint8)
    return np.where(np.isfinite(values), classes, 255).astype(np.uint8)


def sensitivity_engines(engine: ScoringEngine) -> dict[str, ScoringEngine]:
    sensitivity = engine.config.params["sensitivity"]
    variants = {}
    for mean, minimum in sensitivity["mean_min_weights"]:
        params = copy.deepcopy(engine.config.params)
        params["soil_quality"] = {"mean_weight": mean, "min_weight": minimum}
        variants[f"mean_min_{mean:g}_{minimum:g}"] = params
    for floor in sensitivity["fertility_floors"]:
        params = copy.deepcopy(engine.config.params)
        params["viability"]["fertility_floor"] = floor
        variants[f"fertility_floor_{floor:g}"] = params
    for weight in sensitivity["drought_weights"]:
        params = copy.deepcopy(engine.config.params)
        params["drought"]["weight"] = weight
        variants[f"drought_weight_{weight:g}"] = params
    return {name: ScoringEngine(replace(engine.config, params=validate_params(params)), engine.calibration) for name, params in variants.items()}


def class_changes(baseline, changed, inside, cutoffs) -> dict:
    a, b = (score_classes(value, cutoffs) for value in (baseline, changed))
    if a.shape != b.shape or a.shape != np.asarray(inside).shape:
        raise ValueError("class comparison and country mask shapes must agree")
    both = inside & (a != 255) & (b != 255)
    class_count = len(cutoffs) + 1
    matrix = np.bincount(a[both].astype(np.int64) * class_count + b[both], minlength=class_count ** 2).reshape(class_count, class_count)
    return {"compared_cells": int(np.count_nonzero(both)), "changed_class_cells": int(np.count_nonzero(both & (a != b))),
            "became_missing_cells": int(np.count_nonzero(inside & (a != 255) & (b == 255))),
            "became_valid_cells": int(np.count_nonzero(inside & (a == 255) & (b != 255))), "transition_counts": matrix.tolist()}


def _merge_changes(total: dict, increment: dict) -> None:
    for key in ("compared_cells", "changed_class_cells", "became_missing_cells", "became_valid_cells"):
        total[key] = total.get(key, 0) + increment[key]
    total["transition_counts"] = (np.asarray(total.get("transition_counts", np.zeros_like(increment["transition_counts"])), dtype=np.int64) + increment["transition_counts"]).tolist()


def crop_correlations(viability: dict[str, np.ndarray], inside: np.ndarray) -> list[dict]:
    viability = {name: decode(values) for name, values in viability.items()}
    names = list(viability)
    masks = {name: inside & np.isfinite(viability[name]) for name in names}
    common_masks = bool(names) and all(np.array_equal(masks[name], masks[names[0]]) for name in names)
    ranks = {}
    if common_masks:
        for name in names:
            values = viability[name][masks[name]]
            if values.size >= 2 and values.min() != values.max():
                ranks[name] = rankdata(values)
    results = []
    for a, b in combinations(names, 2):
        valid = masks[a] & masks[b]
        count = int(np.count_nonzero(valid))
        correlation = None
        if count >= 2:
            if common_masks:
                if a in ranks and b in ranks:
                    correlation = float(np.corrcoef(ranks[a], ranks[b])[0, 1])
            else:
                va, vb = viability[a][valid], viability[b][valid]
                if va.min() != va.max() and vb.min() != vb.max():
                    correlation = float(spearmanr(va, vb).statistic)
        results.append({"crops": [a, b], "paired_valid_cells": count, "spearman": correlation,
                        "above_specification_0_95": correlation is not None and correlation > 0.95,
                        "reason": None if correlation is not None else "Insufficient paired cells or constant map"})
    return results


def limiting_shares(values: np.ndarray, lookup: dict[str, int], gates: set[str] | None = None) -> dict:
    values = np.asarray(values, dtype=np.float32)
    valid = np.isfinite(values) & (values != grid.NODATA)
    if not np.all(np.isin(values[valid], tuple(lookup.values()))):
        raise ValueError("limiting shares require legal frozen factor IDs")
    total = int(np.count_nonzero(valid))
    factors = []
    additional_count = 0
    for name, index in lookup.items():
        count = int(np.count_nonzero(values == index))
        group = None if gates is None else "gate" if name in gates else "additional"
        if group == "additional":
            additional_count += count
        factors.append({"id": index, "factor": name, "group": group, "cells": count, "share": count / total if total else None})
    return {"valid_limiting_cells": total, "additional_limiting_cells": additional_count if gates is not None else None,
            "additional_limiting_share": additional_count / total if gates is not None and total else None, "factors": factors}


def edge_continuity(inside, confidence, weight, basic: dict, blended: dict, *, fertility_floor: float) -> dict:
    """Separate spatial variation from the extra effect of confidence blending."""
    inside = np.asarray(inside, dtype=bool)
    confidence, weight = decode(confidence), decode(weight)
    basic = {name: decode(values) for name, values in basic.items()}
    blended = {name: decode(values) for name, values in blended.items()}
    if inside.ndim != 2 or basic.keys() != blended.keys() or any(values.shape != inside.shape for values in (confidence, weight, *basic.values(), *blended.values())):
        raise ValueError("edge diagnostics require matching two-dimensional bands and mask")
    positive = confidence > 0
    pairs = []
    for a, b in (((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
                 ((slice(None, -1), slice(None)), (slice(1, None), slice(None)))):
        cross = inside[a] & inside[b] & (positive[a] != positive[b])
        a_pos = positive[a][cross]
        pairs.append((a, b, cross, a_pos))
    count = sum(int(np.count_nonzero(cross)) for _, _, cross, _ in pairs)
    results = {}
    for name in basic:
        residuals, raw_jumps, basic_jumps, boundary_weights = [], [], [], []
        for a, b, cross, a_pos in pairs:
            ba, bb = basic[name][a][cross], basic[name][b][cross]
            sa, sb = blended[name][a][cross], blended[name][b][cross]
            wa, wb = weight[a][cross], weight[b][cross]
            valid = np.isfinite(ba) & np.isfinite(bb) & np.isfinite(sa) & np.isfinite(sb)
            sign = np.where(a_pos, 1, -1)
            raw_jump = (sa.astype(np.float64) - sb) * sign
            basic_jump = (ba.astype(np.float64) - bb) * sign
            residuals.append((raw_jump - basic_jump)[valid])
            raw_jumps.append(raw_jump[valid])
            basic_jumps.append(basic_jump[valid])
            boundary_weights.append(np.where(a_pos, wa, wb)[valid])
        effect = np.concatenate(residuals) if residuals else np.array([])
        boundary_weight = np.concatenate(boundary_weights) if boundary_weights else np.array([])
        multiplier = 1 if name == "soil_quality" else 1 - fertility_floor
        if np.any(np.abs(effect) > multiplier * boundary_weight + ROUNDING_TOLERANCE):
            raise ValueError(f"{name}: confidence-edge effect exceeds the theoretical blend bound")
        zero_weight = inside & (weight == 0) & np.isfinite(basic[name]) & np.isfinite(blended[name])
        mismatches = int(np.count_nonzero(basic[name][zero_weight].view(np.uint32) != blended[name][zero_weight].view(np.uint32)))
        if mismatches:
            raise ValueError(f"{name}: zero-weight scores do not exactly equal basic scores")
        results[name] = {"valid_boundary_pairs": int(effect.size), "absolute_confidence_only_jump": numeric_summary(np.abs(effect)),
                         "absolute_raw_spatial_jump": numeric_summary(np.abs(np.concatenate(raw_jumps))) if raw_jumps else numeric_summary([]),
                         "absolute_basic_spatial_jump": numeric_summary(np.abs(np.concatenate(basic_jumps))) if basic_jumps else numeric_summary([]),
                         "positive_side_weight": numeric_summary(boundary_weight), "theoretical_bound_multiplier": multiplier,
                         "zero_weight_compared_cells": int(np.count_nonzero(zero_weight)), "zero_weight_bitwise_mismatches": mismatches}
    return {"horizontal_and_vertical_coverage_boundary_pairs": count, "bands": results,
            "method": "Inside-country neighbors crossing supplied positive vs zero EthioSIS confidence. Subtract the basic spatial jump from the blended jump; no artificial interpolation or smoothing.",
            "interpretation": "Exact zero weight agrees with basic and refinement is bounded by its weight. These numerical diagnostics do not replace owner/visual review of the source coverage boundary."}


def _middle_cell(mask: np.ndarray) -> tuple[int, int] | None:
    candidates = np.flatnonzero(mask)
    return tuple(int(value) for value in np.unravel_index(candidates[len(candidates) // 2], mask.shape)) if candidates.size else None


def _nearest_cell(mask: np.ndarray, row: int, col: int) -> tuple[int, int] | None:
    rows, cols = np.nonzero(mask)
    if not rows.size:
        return None
    # Diagnostic selection uses grid-index distance, NOT a road/travel distance.
    # np.nonzero is row-major, so equal-distance ties are deterministic.
    selected = int(np.argmin((rows.astype(np.int64) - row) ** 2 + (cols.astype(np.int64) - col) ** 2))
    return int(rows[selected]), int(cols[selected])


def geographic_checks(source, engine: ScoringEngine, inside, viability: dict, elevation, rain_means: dict, references: dict | None) -> dict:
    indices = stack.band_indices_by_name(source)
    eligible = inside.copy()
    for values in viability.values():
        eligible &= np.isfinite(values)
    locations = {}
    reference_context = {}
    reference_records = references.get("geographic_checks", {}).get("sites", {}) if references else {}
    for source_name, label in (("Addis Ababa", "Central highlands near Addis Ababa"),
                               ("Danakil/Dallol diagnostic site", "Afar/Danakil lowlands")):
        if source_name in reference_records:
            ref = reference_records[source_name]
            row, col = int(ref["common_grid_row"]), int(ref["common_grid_col"])
            if not 0 <= row < inside.shape[0] or not 0 <= col < inside.shape[1]:
                raise ValueError("Phase 1 diagnostic reference is outside the canonical grid")
            nearest = _nearest_cell(eligible, row, col)
            locations[f"Original {source_name} reference"] = ((row, col), "Unmodified Phase 1 reference cell, including any exclusion or source gaps")
            reference_context[f"Original {source_name} reference"] = {"reference_row": row, "reference_col": col, "grid_distance_cells": 0.0}
            if nearest is not None:
                locations[label] = (nearest, f"Nearest cell with all crop scores known to the Phase 1 {source_name} reference, by grid-index distance; no scores were optimized")
                reference_context[label] = {"reference_row": row, "reference_col": col,
                                            "grid_distance_cells": math.hypot(nearest[0] - row, nearest[1] - col),
                                            "distance_label": "Grid-index diagnostic distance, not exact kilometres, road distance or travel time"}
    high = _middle_cell(eligible & np.isfinite(elevation) & (elevation > 3000))
    if high:
        locations["Above 3000 m"] = (high, "Middle eligible cell in row-major order with elevation >3000 m, the specification's diagnostic condition")
    sorghum = next((crop for crop in engine.config.crops if crop.crop == "sorghum"), None)
    if sorghum is not None:
        rule = engine.config.requirements["sorghum"]["rain_season_mm"]
        adequate = np.zeros(inside.shape, dtype=bool)
        for season in sorghum.seasons:
            adequate |= (rain_means[season] >= rule.b) & (rain_means[season] <= rule.c)
        low = _middle_cell(eligible & np.isfinite(elevation) & (elevation < 1200) & adequate)
        if low:
            locations["Below 1200 m with adequate rain"] = (low, "Middle eligible cell in row-major order below 1200 m; at least one allowed season is in the owner-supplied sorghum rainfall plateau")
    sites = {}
    for name, ((row, col), selection) in locations.items():
        values = {band: pixels[0, 0] for band, pixels in read_named(source, indices, stack.CORE_BAND_NAMES, (col, row, 1, 1)).items()}
        report = engine.score_cell(values)
        sites[name] = {"row": row, "col": col, "cell_centre_lon_lat": list(gdal.ApplyGeoTransform(grid.TRANSFORM.to_gdal(), col + 0.5, row + 0.5)),
                       "selection": selection, "resolution_label": "Nominal 1 km cell; not an exact 90 m reading",
                       "elevation_m": None if values["elevation_m"] == grid.NODATA else float(values["elevation_m"]),
                       "landcover": None if values["landcover"] == grid.NODATA else int(values["landcover"]),
                       "reference_context": reference_context.get(name), "scores": report["scores"], "crop_breakdowns": report["crops"]}

    def compare(site_a, crop_a, site_b, crop_b):
        a = sites.get(site_a, {}).get("scores", {}).get(f"viability_{crop_a}")
        b = sites.get(site_b, {}).get("scores", {}).get(f"viability_{crop_b}")
        return a > b if a is not None and b is not None else None

    checks = {"highland_teff_above_afar": compare("Central highlands near Addis Ababa", "teff", "Afar/Danakil lowlands", "teff"),
              "highland_wheat_above_afar": compare("Central highlands near Addis Ababa", "wheat", "Afar/Danakil lowlands", "wheat"),
              "above_3000_barley_above_maize": compare("Above 3000 m", "barley", "Above 3000 m", "maize"),
              "above_3000_barley_above_sorghum": compare("Above 3000 m", "barley", "Above 3000 m", "sorghum"),
              "below_1200_sorghum_above_barley": compare("Below 1200 m with adequate rain", "sorghum", "Below 1200 m with adequate rain", "barley"),
              "below_1200_sorghum_above_wheat": compare("Below 1200 m with adequate rain", "sorghum", "Below 1200 m with adequate rain", "wheat")}
    return {"checks": checks, "sites": sites, "selection_note": "Representative sites are not a guarantee for every cell in an elevation belt. Missing or false expectations remain visible for owner review."}


def review_diagnostics(source, candidate, manifest: dict, *, references: dict | None = None) -> dict:
    """Read verified pixels; the builder may pass the bit-identical assembly.

    That optimized path is allowed only after technical acceptance compares all
    assembly/COG/parent pixels. The diagnostic calculations themselves are the
    same whether reading the COG or the verified band-interleaved intermediate.
    """
    engine = stack_v2.frozen_engine(manifest)
    source_indices, indices = stack.band_indices_by_name(source), stack.band_indices_by_name(candidate)
    inside = collect_band(candidate, indices, "ethiopia_mask") == 1
    stats = {record["name"]: band_statistics(collect_band(candidate, indices, record["name"]), inside, record) for record in manifest["bands"]}
    soil_arrays = {name: decode(collect_band(candidate, indices, name)) for name in ("soil_quality_basic", "soil_quality", "data_weight")}
    viability = {crop.crop: decode(collect_band(candidate, indices, f"viability_{crop.crop}")) for crop in engine.config.crops}
    base_crop = {name: np.full(inside.shape, np.nan, dtype=np.float32) for name in viability}
    missing_inputs = {name: 0 for name in sorted(engine.required_bands - {"ethiopia_mask"})}
    gate_counts = {crop.crop: {"missing": {}, "zero": {}, "chosen_seasons": {season: 0 for season in crop.seasons},
                              "eligible_landcover_cells": 0, "zero_gate_union_cells": 0, "missing_gate_union_cells": 0,
                              "missing_additional_basic_cells": 0, "slope_hard_zero_cells": 0,
                              "unknown_chosen_season_cells": 0} for crop in engine.config.crops}
    variants = sensitivity_engines(engine)
    sensitivity = {name: {"parameters": variant.config.params, "bands": {}, "best_crop_id_changes": 0} for name, variant in variants.items()}
    compared_names = ("soil_quality_basic", "soil_quality_detailed", "soil_quality", "soil_score_overall", *(f"viability_{crop.crop}" for crop in engine.config.crops))
    for window in windows(source.RasterXSize, source.RasterYSize):
        col, row, width, height = window
        values = read_named(source, source_indices, stack.CORE_BAND_NAMES, window)
        mask = values["ethiopia_mask"] == 1
        evaluated = engine.evaluate(values, breakdown=True)
        for name in missing_inputs:
            missing_inputs[name] += int(np.count_nonzero(mask & ~np.isfinite(decode(values[name]))))
        for crop in engine.config.crops:
            detail = evaluated.crop_details[crop.crop]
            gate_scores = [item["score"] for item in detail["factors"].values() if item["group"] == "gate"]
            basic, _, _ = crops.viability(gate_scores, detail["additional_basic"], slope_score=detail["factors"]["slope_deg"]["score"], floor=engine.config.params["viability"]["fertility_floor"])
            allowed = mask & np.isfinite(decode(values["landcover"])) & ~np.isin(values["landcover"], (50, 80))
            base_crop[crop.crop][row:row + height, col:col + width] = np.where(allowed, basic, np.nan)
            for name, condition in {"eligible_landcover_cells": allowed,
                                    "zero_gate_union_cells": allowed & np.any(np.stack(gate_scores) == 0, axis=0),
                                    "missing_gate_union_cells": allowed & np.any(~np.isfinite(np.stack(gate_scores)), axis=0),
                                    "missing_additional_basic_cells": allowed & ~np.isfinite(detail["additional_basic"]),
                                    "slope_hard_zero_cells": allowed & (detail["factors"]["slope_deg"]["score"] == 0),
                                    "unknown_chosen_season_cells": allowed & (detail["chosen_season"] < 0)}.items():
                gate_counts[crop.crop][name] += int(np.count_nonzero(condition))
            for factor, item in detail["factors"].items():
                if item["group"] == "gate":
                    for label, condition in (("missing", ~np.isfinite(item["score"])), ("zero", item["score"] == 0)):
                        gate_counts[crop.crop][label][factor] = gate_counts[crop.crop][label].get(factor, 0) + int(np.count_nonzero(allowed & condition))
            for season_index, season in enumerate(crop.seasons):
                gate_counts[crop.crop]["chosen_seasons"][season] += int(np.count_nonzero(allowed & (detail["chosen_season"] == season_index)))
        for name, variant in variants.items():
            changed = variant.evaluate(values).bands
            for band_name in compared_names:
                increment = class_changes(evaluated.bands[band_name], changed[band_name], mask, engine.config.params["score_classes"])
                _merge_changes(sensitivity[name]["bands"].setdefault(band_name, {}), increment)
            sensitivity[name]["best_crop_id_changes"] += int(np.count_nonzero(mask & (evaluated.bands["best_crop"] != changed["best_crop"])))
    for result in sensitivity.values():
        for band in result["bands"].values():
            band["changed_class_share"] = band["changed_class_cells"] / band["compared_cells"] if band["compared_cells"] else None
    limits = {crop.crop: limiting_shares(collect_band(candidate, indices, f"viability_{crop.crop}_limiting"), engine.crop_lookups[crop.crop], set(engine.crop_rules(crop)[0])) for crop in engine.config.crops}
    confidence = decode(collect_band(source, source_indices, "es_confidence"))
    edge = edge_continuity(inside, confidence, soil_arrays["data_weight"], {"soil_quality": soil_arrays["soil_quality_basic"], **{f"viability_{name}": value for name, value in base_crop.items()}},
                           {"soil_quality": soil_arrays["soil_quality"], **{f"viability_{name}": value for name, value in viability.items()}}, fertility_floor=engine.config.params["viability"]["fertility_floor"])
    elevation = decode(collect_band(source, source_indices, "elevation_m"))
    rain_means = {season: decode(collect_band(source, source_indices, f"{season}_mean_mm")) for season in ("kiremt", "belg")}
    geography = geographic_checks(source, engine, inside, viability, elevation, rain_means, references)
    correlations = crop_correlations(viability, inside)
    weight = soil_arrays["data_weight"]
    cutoff = engine.config.params["blend"]["detailed_label_min_weight"]
    warnings = []
    for crop in engine.config.crops:
        statistic = stats[f"viability_{crop.crop}"]
        if not statistic["valid_inside_cells"]:
            warnings.append({"code": f"no_valid_viability_{crop.crop}", "message": f"{crop.display_name} has no valid inside-country scores"})
        elif statistic["exact_zero_cells"] > statistic["valid_inside_cells"] / 2:
            warnings.append({"code": f"majority_zero_viability_{crop.crop}", "message": f"More than half the valid {crop.display_name} cells have hard-zero viability; review gate-zero counts and ranges"})
    for pair in correlations:
        if pair["above_specification_0_95"]:
            warnings.append({"code": "high_correlation_" + "_".join(pair["crops"]), "message": f"Spearman correlation {pair['spearman']:.6f} exceeds the specification's 0.95 review flag"})
    for name, result in geography["checks"].items():
        if result is not True:
            warnings.append({"code": f"geographic_{name}", "message": "Expected geographic comparison failed or could not be evaluated; see exact site inputs and factor scores"})
    if limits and all(record["additional_limiting_share"] is not None and record["additional_limiting_share"] > 0.5 for record in limits.values()):
        warnings.append({"code": "majority_additional_limiting_all_crops", "message": "Additional factors limit a majority of valid cells for every crop; review gate ranges rather than silently retuning"})
    return {"band_statistics": stats, "missing_active_input_cells": missing_inputs, "crop_gate_diagnostics": gate_counts,
            "crop_correlations": correlations, "crop_limiting_shares": limits,
            "soil_limiting_shares": limiting_shares(collect_band(candidate, indices, "soil_limiting_factor"), engine.soil_lookup),
            "coverage": {"inside_cells": int(np.count_nonzero(inside)), "positive_source_confidence_cells": int(np.count_nonzero(inside & (confidence > 0))),
                         "basic_weight_cells": int(np.count_nonzero(inside & (weight == 0))),
                         "partial_weight_cells": int(np.count_nonzero(inside & (weight > 0) & (weight < cutoff))),
                         "detailed_weight_cells": int(np.count_nonzero(inside & (weight >= cutoff))),
                         "interpretation": "Positive confidence describes interpolated support, not cells with exact measured samples; do not equate this footprint with measured-sample coverage."},
            "edge_continuity": edge, "geographic_checks": geography,
            "sensitivity": {"score_class_cutoffs_percent": engine.config.params["score_classes"], "variants": sensitivity,
                            "class_rule": "Float32 score cutoffs; exact cutoffs enter the higher class; missing comparisons are separately counted.",
                            "coverage_radius": {"status": "skipped_owner_requested", "reason": "No 20 km radius regeneration, expansion or shrink simulation was requested; supplied IDW values remain unchanged."}},
            "warnings": warnings, "review_status": "Owner review required; numerical diagnostics are not agronomic validation"}


def markdown_report(report: dict) -> str:
    """Readable summaries; complete histograms and site factor bars are in JSON."""
    diagnostics, technical = report["diagnostics"], report["technical"]
    lines = ["# Stack v2 scoring build report", "", f"Status: **{report['status']}** — candidate only; not promoted.", "",
             f"Candidate: `{report['raster_file']}`", f"Raster SHA-256: `{report['raster_sha256']}`", f"Manifest SHA-256: `{report['manifest_sha256']}`", "",
             "## Technical acceptance", "", f"- {technical['band_count']} bands; the original {technical['v1_band_count']} bands and their manifest records are unchanged.",
             f"- Grid, finite float32/-9999 storage, legal categorical values, country mask, COG structure and all overview pixels: **passed**.",
             f"- Every scored cell matches the shared evaluator; {technical['scalar_checks']['compared_cells']} representative/random single-cell reports match bit-for-bit.",
             f"- Reproducibility: **{report['reproducibility']['status']}**. Original inputs and stack v1 remain unchanged.", "",
             "## Approved configuration", "", f"- CV cutoffs: `{report['calibration']['cv_low']:.9g}` / `{report['calibration']['cv_high']:.9g}` from selected-season Ethiopian cells.",
             "- Reliability is the approved Normal approximation to rainfall reaching a configured fraction of seasonal mean, not a crop water-need probability.",
             "- Phosphorus stays mg/100kg; EthioSIS stays mg/kg; nitrogen is sg_n only. Slope is the accepted mean-degree band.",
             "- Active thresholds, units, aliases, sources, tuning and original crop IDs are frozen in the manifest. No scientific thresholds are embedded in scoring code.", "",
             "## Skipped items", ""]
    lines += [f"- {row.get('crop', 'Soil')}/{row.get('factor', row.get('kind'))}: {row['reason']}" for row in report["skipped"]] or ["- None."]
    lines += ["- Coverage-radius sensitivity: skipped at the owner's request; no interpolation was regenerated.", "", "## Review flags", ""]
    lines += [f"- **{row['code']}**: {row['message']}" for row in diagnostics["warnings"]] or ["- No automatic flags; owner review is still required."]
    lines += ["", "## Band distributions", "", "All counts/ranges are observed values. Complete per-band histograms (including raw inputs and categorical IDs) are in `score_build_report.json`.", "",
              "| Band | Valid inside | Missing inside | Min | Median | Max | Exactly 0 | Exactly 1 |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    def number(value):
        return "null" if value is None else f"{value:.6g}"
    for name, stats in diagnostics["band_statistics"].items():
        lines.append(f"| {name} | {stats['valid_inside_cells']:,} | {stats['missing_inside_cells']:,} | {number(stats['min'])} | {number(stats['median'])} | {number(stats['max'])} | {stats['exact_zero_cells']:,} | {stats['exact_one_cells']:,} |")
    lines += ["", "### Preserved raw-input sanity outliers", "", "These inherited rough ranges are diagnostics, not scoring thresholds or clipping limits.", ""]
    outliers = []
    for name in stack.CORE_BAND_NAMES:
        stats = diagnostics["band_statistics"][name]
        below, above = stats["below_sanity_range_cells"], stats["above_sanity_range_cells"]
        if below or above:
            outliers.append(f"- `{name}`: {below or 0:,} below / {above or 0:,} above the inherited sanity range; values unchanged.")
    lines += outliers or ["- None detected within supplied sanity bounds; unsupplied bounds remain null."]
    lines += ["", "### Missing active inputs", "", "Counts inside Ethiopia, including excluded land cover; zero counts omitted here but present in JSON.", ""]
    lines += [f"- `{name}`: {count:,}" for name, count in diagnostics["missing_active_input_cells"].items() if count] or ["- No missing active inputs."]
    lines += ["", "## Score histograms", "", "Counts in ten equal-width score bins from 0 to 1; the final bin includes 1. Constant/empty bands are explicitly labelled.", ""]
    for name, stats in diagnostics["band_statistics"].items():
        if name in stack_v2.SOIL_BAND_NAMES or name.startswith("viability_") or name in stack_v2.BEST_CROP_BAND_NAMES:
            if stats["histogram"]["kind"] != "categorical":
                lines.append(f"- `{name}`: {stats['histogram']['kind']} — `{stats['histogram']['counts']}`")
    lines += ["", "## Crop separation and limiting factors", "", "| Crop pair | Paired cells | Spearman | >0.95 flag |", "|---|---:|---:|---|"]
    for pair in diagnostics["crop_correlations"]:
        lines.append(f"| {' / '.join(pair['crops'])} | {pair['paired_valid_cells']:,} | {number(pair['spearman'])} | {pair['above_specification_0_95']} |")
    lines += ["", "| Crop | Valid limiting cells | Additional limiting share |", "|---|---:|---:|"]
    for crop, shares in diagnostics["crop_limiting_shares"].items():
        lines.append(f"| {crop} | {shares['valid_limiting_cells']:,} | {number(shares['additional_limiting_share'])} |")
    lines += ["", "### Gate-zero and missing-data diagnosis", "", "Counts are over nonexcluded, known land-cover cells. Gate zeros, missing inputs and slope zeros can overlap; unknown active inputs still take precedence over a hard zero.", "",
              "| Crop | Eligible land cover | Any gate zero | Missing gate | Missing basic additional | Hard-zero slope |", "|---|---:|---:|---:|---:|---:|"]
    for crop, counts in diagnostics["crop_gate_diagnostics"].items():
        lines.append(f"| {crop} | {counts['eligible_landcover_cells']:,} | {counts['zero_gate_union_cells']:,} | {counts['missing_gate_union_cells']:,} | {counts['missing_additional_basic_cells']:,} | {counts['slope_hard_zero_cells']:,} |")
    lines += ["", "| Crop | Gate factor | Exactly zero | Missing |", "|---|---|---:|---:|"]
    for crop, counts in diagnostics["crop_gate_diagnostics"].items():
        for factor, count in counts["zero"].items():
            lines.append(f"| {crop} | {factor} | {count:,} | {counts['missing'].get(factor, 0):,} |")
    lines += ["", "Per-factor limiting shares and chosen-season counts are in the JSON report.", "", "## Geographic checks", "",
              "Capital cells are not silently treated as cropland. Nearby reference replacements use only input availability and disclosed grid-distance selection, not optimized crop rankings.", "",
              "| Site | Cell (row, col) | Elevation m | Teff | Wheat | Maize | Sorghum | Barley | Chickpea |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name, site in diagnostics["geographic_checks"]["sites"].items():
        scores = " | ".join(number(site["scores"].get(f"viability_{crop}")) for crop in ("teff", "wheat", "maize", "sorghum", "barley", "chickpea"))
        lines.append(f"| {name} | ({site['row']}, {site['col']}) | {number(site['elevation_m'])} | {scores} |")
    lines += [""] + [f"- {name}: **{result}**" for name, result in diagnostics["geographic_checks"]["checks"].items()]
    lines += ["", "## EthioSIS edge continuity", "", "The confidence-only effect below subtracts ordinary basic-tier spatial variation from the raw neighboring score difference. It is bounded by blend weight (and by 1-floor for crop modifiers). Zero-weight cells match basic scores exactly; visual/scientific boundary review remains required.", "",
              "| Band | Valid boundary pairs | Median absolute confidence effect | Max absolute confidence effect | Max positive-side weight |", "|---|---:|---:|---:|---:|"]
    for name, result in diagnostics["edge_continuity"]["bands"].items():
        lines.append(f"| {name} | {result['valid_boundary_pairs']:,} | {number(result['absolute_confidence_only_jump']['median'])} | {number(result['absolute_confidence_only_jump']['max'])} | {number(result['positive_side_weight']['max'])} |")
    lines += ["", "## Approved sensitivity variants", "", f"Score-class cutoffs (%): `{diagnostics['sensitivity']['score_class_cutoffs_percent']}`. Missing-status changes are counted separately in JSON.", "",
              "| Variant | Band | Compared cells | Changed class | Changed share |", "|---|---|---:|---:|---:|"]
    for name, result in diagnostics["sensitivity"]["variants"].items():
        for band, counts in result["bands"].items():
            lines.append(f"| {name} | {band} | {counts['compared_cells']:,} | {counts['changed_class_cells']:,} | {number(counts['changed_class_share'])} |")
    lines += ["", "## Review and next step", "", "Review the gate-zero distributions, geographic expectations, crop correlations, limiting shares and sensitivity before approving this exact candidate. No Django, API, UI, PostGIS or PMTiles work is included.", ""]
    return "\n".join(lines)
