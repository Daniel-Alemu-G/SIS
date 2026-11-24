"""Unit-free factor functions; all scientific parameters are supplied by callers."""

from __future__ import annotations

import math

import numpy as np
from scipy.special import ndtr


def _interval(a: float, b: float) -> None:
    if not all(math.isfinite(v) for v in (a, b)) or not a < b:
        raise ValueError("breakpoints must be finite with a < b")


def _values(x) -> np.ndarray:
    values = np.asarray(x, dtype=np.float32)
    return np.where(np.isfinite(values), values, np.nan).astype(np.float32)


def ramp_up(x, a: float, b: float) -> np.ndarray:
    """Zero at/below a, one at/above b, linear in between; missing stays NaN."""
    _interval(a, b)
    values = _values(x).astype(np.float64)
    return np.clip((values - a) / (b - a), 0, 1).astype(np.float32)


def ramp_down(x, a: float, b: float) -> np.ndarray:
    return (1 - ramp_up(x, a, b)).astype(np.float32)


def trapezoid(x, a: float, b: float, c: float, d: float) -> np.ndarray:
    _interval(a, b)
    _interval(c, d)
    if b > c:
        raise ValueError("trapezoid breakpoints must satisfy a < b <= c < d")
    values = _values(x).astype(np.float64)
    up = np.clip((values - a) / (b - a), 0, 1)
    down = np.clip((d - values) / (d - c), 0, 1)
    return np.minimum(up, down).astype(np.float32)


def rainfall_reliability_normal(cv, *, threshold_ratio: float) -> np.ndarray:
    """Owner-approved P(rain >= ratio * mean), in percent, using a Normal model.

    This is NOT P(rain >= a crop's rainfall requirement). It uses only seasonal
    CV and a configured fraction of mean; no monthly data is regenerated.
    Exact zero CV uses the owner's deterministic convention. Negative or
    non-finite CV is missing, never converted into reliable rainfall.
    """
    if not math.isfinite(threshold_ratio) or threshold_ratio < 0:
        raise ValueError("threshold_ratio must be finite and nonnegative")
    values = np.asarray(cv, dtype=np.float64)
    z = np.full(values.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(values) & (values > 0)
    np.divide(1 - threshold_ratio, values, out=z, where=valid)
    result = ndtr(z) * 100
    result = np.where(values == 0, 100 if threshold_ratio <= 1 else 0, result)
    return np.where(np.isfinite(values) & (values >= 0), result, np.nan).astype(np.float32)
