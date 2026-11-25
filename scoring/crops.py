"""Crop gates, season selection, modifiers and deterministic best-crop logic."""

from __future__ import annotations

import math

import numpy as np

from scoring.soil import score_stack


def choose_season(rainfall_scores, rainfall_means) -> np.ndarray:
    """Best rainfall score, then higher mean, then configured season order.

    All allowed-season means are needed for the comparison. A missing allowed
    season propagates rather than being silently dropped.
    """
    scores = score_stack(rainfall_scores)
    means = np.stack([np.asarray(v, dtype=np.float32) for v in rainfall_means], axis=0)
    if means.shape != scores.shape:
        raise ValueError("rainfall score/mean shapes differ")
    chosen = np.zeros(scores.shape[1:], dtype=np.int8)
    best_score = scores[0]
    best_mean = means[0]
    for index in range(1, scores.shape[0]):
        better = (scores[index] > best_score) | ((scores[index] == best_score) & (means[index] > best_mean))
        chosen = np.where(better, index, chosen).astype(np.int8)
        best_score = np.where(better, scores[index], best_score)
        best_mean = np.where(better, means[index], best_mean)
    valid = np.all(np.isfinite(scores) & np.isfinite(means) & (means >= 0), axis=0)
    return np.where(valid, chosen, -1).astype(np.int8)


def selected_values(arrays, chosen) -> np.ndarray:
    values = np.stack([np.asarray(array, dtype=np.float32) for array in arrays], axis=0)
    if values.shape[1:] != np.asarray(chosen).shape:
        raise ValueError("season choice and values have different shapes")
    result = np.take_along_axis(values, np.maximum(chosen, 0)[None, ...], axis=0)[0]
    return np.where(chosen >= 0, result, np.nan).astype(np.float32)


def hybrid_gate(gate_scores) -> np.ndarray:
    values = score_stack(gate_scores).astype(np.float64)
    any_missing = np.any(np.isnan(values), axis=0)
    any_zero = np.any(values == 0, axis=0)
    # Positive values use logs to avoid product underflow. Replacing zero by
    # one is a numerical device only; any zero is restored as a hard zero.
    log_values = np.log(np.where(values > 0, values, 1))
    geometric = np.exp(np.mean(log_values, axis=0))
    minimum = np.min(values, axis=0)
    gate = np.sqrt(minimum * geometric)
    return np.where(any_missing, np.nan, np.where(any_zero, 0, gate)).astype(np.float32)


def viability(gate_scores, additional_score, *, slope_score, floor: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not math.isfinite(floor) or not 0 <= floor <= 1:
        raise ValueError("fertility floor must be in [0, 1]")
    gate = hybrid_gate(gate_scores)
    additional, slope = (np.asarray(v, dtype=np.float32) for v in (additional_score, slope_score))
    if additional.shape != gate.shape or slope.shape != gate.shape:
        raise ValueError("gate/additional/slope shapes differ")
    score_stack([additional, slope])
    modifier = (floor + (1 - floor) * additional.astype(np.float64)).astype(np.float32)
    out = (gate.astype(np.float64) * modifier).astype(np.float32)
    # Missing active gate/basic input takes precedence over hard zeros. This
    # follows the no-imputation rule: unknown inputs never become fake zero.
    known = np.isfinite(gate) & np.isfinite(additional) & np.isfinite(slope)
    out = np.where(known, np.where(slope == 0, 0, out), np.nan).astype(np.float32)
    return out, gate, modifier


def best_crop(scores: dict[str, np.ndarray], crop_indices: dict[str, int], *, minimum_score: float) -> dict[str, np.ndarray]:
    """Ignore missing crops, keep original 1-based IDs, and use stable tie order.

    Internal logical nodata is 255. The physical stack writer converts it to
    the common -9999 encoding. A sole valid crop has no known runner-up, so its
    margin is NaN rather than an invented comparison with zero.
    """
    if not scores:
        raise ValueError("best_crop requires at least one computed crop")
    if not math.isfinite(minimum_score) or not 0 <= minimum_score <= 1:
        raise ValueError("best-crop minimum score must be in [0, 1]")
    names = sorted(scores, key=lambda name: crop_indices[name])
    if len(set(crop_indices[name] for name in names)) != len(names) or any(not 1 <= crop_indices[name] < 255 for name in names):
        raise ValueError("crop indices must be unique in 1..254")
    values = score_stack([scores[name] for name in names])
    finite = np.isfinite(values)
    valid_count = np.sum(finite, axis=0)
    ranked = np.where(finite, values, -np.inf)
    winner = np.argmax(ranked, axis=0)
    best = np.take_along_axis(ranked, winner[None, ...], axis=0)[0]
    ids = np.asarray([crop_indices[name] for name in names], dtype=np.uint8)[winner]
    ids = np.where(valid_count == 0, 255, np.where(best < minimum_score, 0, ids)).astype(np.uint8)
    margin = np.full(best.shape, np.nan, dtype=np.float32)
    if len(names) > 1:
        runner_up = np.sort(ranked, axis=0)[-2]
        np.subtract(best, runner_up, out=margin, where=valid_count >= 2)
    return {"best_crop": ids, "best_crop_score": np.where(valid_count, best, np.nan).astype(np.float32),
            "best_crop_margin": margin}
