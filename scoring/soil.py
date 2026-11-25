"""Soil aggregation, confidence blending and data-derived drought calibration."""

from __future__ import annotations

import math

import numpy as np

from scoring.functions import ramp_up


def score_stack(scores) -> np.ndarray:
    if not scores:
        raise ValueError("at least one factor score is required")
    arrays = [np.asarray(score, dtype=np.float32) for score in scores]
    if any(array.shape != arrays[0].shape for array in arrays):
        raise ValueError("all factor scores must have the same shape")
    result = np.stack(arrays, axis=0)
    finite = np.isfinite(result)
    if np.any(np.isinf(result)) or np.any(finite & ((result < 0) | (result > 1))):
        raise ValueError("factor scores must be in [0, 1] or NaN")
    return result


def soil_quality(basic_scores, extra_es_scores=None, *, mean_weight: float, min_weight: float) -> np.ndarray:
    """Basic missing factors propagate; only extra EthioSIS missingness is ignored."""
    if not all(math.isfinite(w) and 0 <= w <= 1 for w in (mean_weight, min_weight)) or not math.isclose(mean_weight + min_weight, 1, abs_tol=1e-12):
        raise ValueError("mean/min weights must be nonnegative and sum to 1")
    basic = score_stack(basic_scores)
    basic_missing = np.any(np.isnan(basic), axis=0)
    combined = basic
    no_es = np.zeros(basic.shape[1:], dtype=bool)
    if extra_es_scores:
        es = score_stack(extra_es_scores)
        if es.shape[1:] != basic.shape[1:]:
            raise ValueError("basic and EthioSIS factors must have the same shape")
        no_es = np.all(np.isnan(es), axis=0)
        combined = np.concatenate((basic, es), axis=0)
    present = np.isfinite(combined)
    count = np.sum(present, axis=0)
    total = np.sum(np.where(present, combined, 0), axis=0, dtype=np.float64)
    mean = np.full(count.shape, np.nan, dtype=np.float64)
    np.divide(total, count, out=mean, where=count > 0)
    minimum = np.min(np.where(present, combined, np.inf), axis=0)
    result = mean_weight * mean + min_weight * minimum
    return np.where(basic_missing | no_es, np.nan, result).astype(np.float32)


def data_weight(confidence, extra_es_scores) -> np.ndarray:
    """Confidence times active factor completeness; missing confidence means zero."""
    values = np.asarray(confidence, dtype=np.float32)
    if np.any(np.isfinite(values) & ((values < 0) | (values > 1))) or np.any(np.isinf(values)):
        raise ValueError("confidence must be in [0, 1] or NaN")
    clean = np.where(np.isnan(values), 0, values)
    if not extra_es_scores:
        return np.zeros(values.shape, dtype=np.float32)
    scores = score_stack(extra_es_scores)
    if scores.shape[1:] != values.shape:
        raise ValueError("confidence and factor shapes differ")
    fraction = np.sum(np.isfinite(scores), axis=0) / scores.shape[0]
    return (clean.astype(np.float64) * fraction).astype(np.float32)


def blend(basic, detailed, weight) -> np.ndarray:
    basic, detailed, weight = (np.asarray(v, dtype=np.float32) for v in (basic, detailed, weight))
    if basic.shape != detailed.shape or basic.shape != weight.shape:
        raise ValueError("blend inputs must have identical shapes")
    if np.any(np.isfinite(weight) & ((weight < 0) | (weight > 1))) or np.any(np.isinf(weight)):
        raise ValueError("blend weight must be in [0, 1] or NaN")
    use_basic = np.isnan(detailed) | np.isnan(weight) | (weight <= 0)
    result = (1 - weight.astype(np.float64)) * basic + weight.astype(np.float64) * detailed
    return np.where(use_basic, basic, result).astype(np.float32)


def limiting_factor(scores: dict, lookup: dict[str, int], *, valid=None) -> np.ndarray:
    """First lowest factor in the frozen lookup order wins an exact tie."""
    names = list(scores)
    values = score_stack(list(scores.values()))
    if any(name not in lookup for name in names):
        raise ValueError("factor lookup is incomplete")
    present = np.isfinite(values)
    selected = np.argmin(np.where(present, values, np.inf), axis=0)
    ids = np.asarray([lookup[name] for name in names], dtype=np.float32)[selected]
    usable = np.any(present, axis=0)
    if valid is not None:
        usable &= np.asarray(valid, dtype=bool)
    return np.where(usable, ids, np.nan).astype(np.float32)


def main_rainy_season(kiremt_mean, belg_mean, kiremt_cv, belg_cv) -> tuple[np.ndarray, np.ndarray]:
    """Select the larger rainy-season mean; exact ties choose kiremt.

    Both means are required to make this comparison. Missing comparison data
    is not silently replaced by whichever season happened to be available.
    Codes are 1=kiremt, 2=belg, 0=unknown; bega is never considered.
    """
    k, b, kc, bc = [np.asarray(v, dtype=np.float32) for v in (kiremt_mean, belg_mean, kiremt_cv, belg_cv)]
    if any(v.shape != k.shape for v in (b, kc, bc)):
        raise ValueError("rainfall inputs must have the same shape")
    valid_mean = np.isfinite(k) & np.isfinite(b) & (k >= 0) & (b >= 0)
    choose_k = k >= b
    season = np.where(valid_mean, np.where(choose_k, 1, 2), 0).astype(np.uint8)
    cv = np.where(choose_k, kc, bc)
    cv = np.where(valid_mean & np.isfinite(cv) & (cv >= 0), cv, np.nan).astype(np.float32)
    return season, cv


def calibrate_drought(kiremt_mean, belg_mean, kiremt_cv, belg_cv, inside, *, low_percentile: float, high_percentile: float, method: str) -> dict:
    _, selected_cv = main_rainy_season(kiremt_mean, belg_mean, kiremt_cv, belg_cv)
    inside = np.asarray(inside, dtype=bool)
    if inside.shape != selected_cv.shape:
        raise ValueError("calibration mask shape differs from rainfall")
    sample = selected_cv[inside & np.isfinite(selected_cv)].astype(np.float64)
    return calibrate_selected_cv(sample, missing_cells=int(np.count_nonzero(inside & ~np.isfinite(selected_cv))),
                                 low_percentile=low_percentile, high_percentile=high_percentile, method=method)


def calibrate_selected_cv(sample, *, missing_cells: int, low_percentile: float, high_percentile: float, method: str) -> dict:
    """Calibrate collected window samples identically to the whole-array helper."""
    if not 0 <= low_percentile < high_percentile <= 100 or method != "linear":
        raise ValueError("CV percentiles must be ordered in 0..100 with explicit linear interpolation")
    sample = np.asarray(sample, dtype=np.float64).ravel()
    if not np.all(np.isfinite(sample) & (sample >= 0)) or type(missing_cells) is not int or missing_cells < 0:
        raise ValueError("selected CV samples must be finite/nonnegative and missing count valid")
    if not sample.size:
        raise ValueError("no valid selected-season CV cells inside Ethiopia")
    low, high = np.percentile(sample, [low_percentile, high_percentile], method=method)
    if not math.isfinite(low) or not math.isfinite(high) or low >= high:
        raise ValueError("CV calibration produced equal/invalid cut-offs")
    return {"cv_low": float(low), "cv_high": float(high), "low_percentile": float(low_percentile),
            "high_percentile": float(high_percentile), "percentile_method": method,
            "valid_calibration_cells": int(sample.size),
            "missing_calibration_cells": missing_cells,
            "selection": "CV of the larger kiremt/belg mean; exact mean ties choose kiremt; bega excluded"}


def drought_score(quality, selected_cv, *, cv_low: float, cv_high: float, weight: float) -> tuple[np.ndarray, np.ndarray]:
    if not math.isfinite(weight) or not 0 <= weight <= 1:
        raise ValueError("drought weight must be in [0, 1]")
    risk = ramp_up(selected_cv, cv_low, cv_high)
    result = np.asarray(quality, dtype=np.float64) * (1 - weight * risk.astype(np.float64))
    return risk, result.astype(np.float32)
