"""Strict, manifest-backed loading of owner thresholds, crops and parameters.

Free-text legacy skip/unit notes are not a machine-readable policy. The
approved YAML explicitly names active/skipped factors and confirmed units;
the original CSVs stay untouched. PROVISIONAL/TBD sources remain unscored.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import yaml

from scoring import functions

SOIL_BINDINGS = {
    "ph": "sg_phh2o", "organic_carbon": "sg_soc", "nitrogen": "sg_n",
    "phosphorus": "sg_p", "cec": "sg_cec", "coarse_fragments": "sg_cfvo",
    "slope": "slope_mean_deg", "potassium": "es_k", "magnesium": "es_mg",
    "sulfur": "es_s", "zinc": "es_zn", "copper": "es_cu", "iron": "es_fe",
    "calcium": "es_ca", "manganese": "es_mn", "molybdenum": "es_mo", "sodium": "es_na",
}
GATE_ORDER = ("elevation_m", "rain_season_mm", "temperature_season_c", "rain_reliability", "ph", "clay_pct")
ADDITIONAL_BASIC = ("slope_deg", "nitrogen", "phosphorus", "organic_carbon", "cec")
SEASONS = ("kiremt", "belg")
SHAPES = {"ramp_up", "ramp_down", "trapezoid"}
SOIL_COLUMNS = {"factor", "display_name", "band", "detailed_band", "shape", "a", "b", "c", "d", "unit", "source", "low_text", "high_text", "notes"}
CROP_COLUMNS = {"crop", "factor", "shape", "a", "b", "c", "d", "unit", "source", "notes"}


class ConfigurationError(ValueError):
    """An ambiguous, malformed or incompatible configuration cannot be used."""


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def canonical_unit(unit: str) -> str:
    # Spelling aliases only. This performs no rescaling or unit conversion.
    return {"°C": "degrees Celsius", "C": "degrees Celsius", "degrees": "degrees",
            "fraction": "fraction", "%": "%"}.get(unit.strip(), unit.strip())


def _number(value, name: str) -> float:
    if isinstance(value, bool):
        raise ConfigurationError(f"{name}: expected a finite number, not boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name}: expected a finite number") from exc
    if not math.isfinite(result):
        raise ConfigurationError(f"{name}: expected a finite number")
    return result


def _bounded(value, name: str, low: float, high: float, *, include_low=True) -> float:
    number = _number(value, name)
    if not ((number >= low if include_low else number > low) and number <= high):
        raise ConfigurationError(f"{name}: outside permitted range")
    return number


def _rows(path: Path, columns: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or set(header) != columns:
            raise ConfigurationError(f"{path}: missing, duplicate or unknown CSV headers")
        result = []
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ConfigurationError(f"{path}:{line}: incorrect CSV field count")
            result.append({key: value.strip() for key, value in row.items()})
    if not result:
        raise ConfigurationError(f"{path}: empty CSV")
    return result


@dataclass(frozen=True)
class Factor:
    factor: str
    shape: str
    a: float
    b: float
    c: float | None
    d: float | None
    unit: str
    source: str
    band: str = ""
    detailed_band: str = ""
    display_name: str = ""
    low_text: str = ""
    high_text: str = ""
    notes: str = ""
    csv_factor: str = ""

    @property
    def input_band(self) -> str:
        return self.band or self.detailed_band

    def score(self, values) -> np.ndarray:
        if self.shape == "trapezoid":
            return functions.trapezoid(values, self.a, self.b, self.c, self.d)
        return getattr(functions, self.shape)(values, self.a, self.b)

    def metadata(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Crop:
    crop: str
    display_name: str
    seasons: tuple[str, ...]
    index: int


@dataclass(frozen=True)
class ScoringConfig:
    soil_factors: tuple[Factor, ...]
    crops: tuple[Crop, ...]
    requirements: dict[str, dict[str, Factor]]
    params: dict
    band_units: dict[str, str]
    skipped: tuple[dict, ...] = ()
    inputs: dict[str, dict] = field(default_factory=dict)

    @property
    def basic_factors(self) -> tuple[Factor, ...]:
        return tuple(factor for factor in self.soil_factors if factor.band)

    @property
    def es_factors(self) -> tuple[Factor, ...]:
        return tuple(factor for factor in self.soil_factors if factor.detailed_band)

    def metadata(self) -> dict:
        return {"params": self.params, "soil_factors": [f.metadata() for f in self.soil_factors],
                "crops": [asdict(crop) for crop in self.crops],
                "requirements": {crop: {name: factor.metadata() for name, factor in rows.items()}
                                 for crop, rows in self.requirements.items()},
                "skipped": list(self.skipped), "inputs": self.inputs}


def validate_params(params: dict) -> dict:
    if not isinstance(params, dict):
        raise ConfigurationError("params.yaml must be a mapping")
    required = {"soil_quality", "drought", "blend", "viability", "score_classes", "rainfall_reliability",
                "crop_factor_aliases", "active_soil_factors", "skipped_soil_factors", "confirmed_units", "approval", "sensitivity"}
    if set(params) != required:
        raise ConfigurationError("params.yaml has missing or unknown sections")
    for key in required - {"score_classes", "active_soil_factors"}:
        if not isinstance(params[key], dict):
            raise ConfigurationError(f"{key}: expected a mapping")
    schemas = {"soil_quality": {"mean_weight", "min_weight"},
               "drought": {"weight", "cv_low_percentile", "cv_high_percentile", "percentile_method"},
               "blend": {"detailed_label_min_weight"},
               "viability": {"fertility_floor", "best_crop_min_score", "best_crop_tie_margin"},
               "rainfall_reliability": {"method", "threshold_ratio"},
               "sensitivity": {"mean_min_weights", "fertility_floors", "drought_weights", "coverage_radius"}}
    for section, keys in schemas.items():
        if set(params[section]) != keys:
            raise ConfigurationError(f"{section}: missing or unknown parameters")
    weights = params["soil_quality"]
    for key in ("mean_weight", "min_weight"):
        _bounded(weights[key], key, 0, 1)
    if not math.isclose(sum(float(weights[k]) for k in weights), 1, abs_tol=1e-12):
        raise ConfigurationError("soil quality weights must sum to 1")
    drought = params["drought"]
    _bounded(drought["weight"], "drought.weight", 0, 1)
    low = _bounded(drought["cv_low_percentile"], "cv_low_percentile", 0, 100)
    high = _bounded(drought["cv_high_percentile"], "cv_high_percentile", 0, 100)
    if low >= high or drought["percentile_method"] != "linear":
        raise ConfigurationError("CV percentiles must be ordered and use explicit linear interpolation")
    _bounded(params["blend"]["detailed_label_min_weight"], "detailed_label_min_weight", 0, 1, include_low=False)
    for key, value in params["viability"].items():
        _bounded(value, key, 0, 1)
    classes = params["score_classes"]
    if not isinstance(classes, list) or not classes:
        raise ConfigurationError("score_classes must be a nonempty ordered list")
    levels = [_bounded(v, "score class", 0, 100, include_low=False) for v in classes]
    if levels[-1] >= 100 or any(a >= b for a, b in zip(levels, levels[1:])):
        raise ConfigurationError("score_classes must be strictly increasing below 100")
    reliability = params["rainfall_reliability"]
    if reliability["method"] not in {"normal_fraction_of_mean", "cv"}:
        raise ConfigurationError("unknown rainfall reliability method")
    _number(reliability["threshold_ratio"], "threshold_ratio")
    if float(reliability["threshold_ratio"]) < 0:
        raise ConfigurationError("threshold_ratio must be nonnegative")
    aliases = params["crop_factor_aliases"]
    if any(key not in {"temperature", "temprature"} or value != "temperature_season_c" for key, value in aliases.items()):
        raise ConfigurationError("unknown or ambiguous crop factor alias")
    active = params["active_soil_factors"]
    skipped = params["skipped_soil_factors"]
    if not isinstance(active, list) or any(not isinstance(v, str) for v in active) or len(active) != len(set(active)):
        raise ConfigurationError("active_soil_factors must be a unique list")
    if set(active) & set(skipped) or (set(active) | set(skipped)) - set(SOIL_BINDINGS):
        raise ConfigurationError("unknown or conflicting soil activation policy")
    if any(not isinstance(reason, str) or not reason.strip() for reason in skipped.values()):
        raise ConfigurationError("every skipped factor needs a reason")
    if set(params["confirmed_units"]) - set(SOIL_BINDINGS):
        raise ConfigurationError("confirmed_units contains an unknown factor")
    if any(not isinstance(unit, str) or not unit.strip() for unit in params["confirmed_units"].values()):
        raise ConfigurationError("confirmed_units values must be unit strings")
    approval = params["approval"]
    if approval.get("approved") is not True or not isinstance(approval.get("authority"), str) or not approval["authority"].strip():
        raise ConfigurationError("owner approval and authority are required")
    if not isinstance(approval.get("chickpea_rainfall_proxy"), str) or not approval["chickpea_rainfall_proxy"].strip():
        raise ConfigurationError("chickpea rainfall proxy must be documented")
    sensitivity = params["sensitivity"]
    if sensitivity["coverage_radius"] != "skip_owner_requested":
        raise ConfigurationError("coverage-radius regeneration is outside this scoring build")
    for section in ("mean_min_weights", "fertility_floors", "drought_weights"):
        if not isinstance(sensitivity[section], list):
            raise ConfigurationError(f"sensitivity.{section} must be a list")
    for pair in sensitivity["mean_min_weights"]:
        if not isinstance(pair, list) or len(pair) != 2:
            raise ConfigurationError("sensitivity mean/min variants must be pairs")
        for number in pair:
            _bounded(number, "sensitivity weight", 0, 1)
        if not math.isclose(sum(pair), 1, abs_tol=1e-12):
            raise ConfigurationError("sensitivity weights must sum to 1")
    for section in ("fertility_floors", "drought_weights"):
        for value in sensitivity[section]:
            _bounded(value, f"sensitivity.{section}", 0, 1)
    return params


def _factor(row: dict, factor: str, *, band="", detailed_band="", template: Factor | None = None) -> Factor:
    shape = row["shape"]
    if shape not in SHAPES:
        raise ConfigurationError(f"{factor}: unknown shape {shape!r}")
    a, b = (_number(row[key], f"{factor}.{key}") for key in ("a", "b"))
    c = d = None
    if shape == "trapezoid":
        c, d = (_number(row[key], f"{factor}.{key}") for key in ("c", "d"))
        if not a < b <= c < d:
            raise ConfigurationError(f"{factor}: require a < b <= c < d")
    elif not a < b or row["c"] or row["d"]:
        raise ConfigurationError(f"{factor}: ramps require a < b and empty c/d")
    name = row.get("display_name") or (template.display_name if template else factor.replace("_", " ").capitalize())
    return Factor(factor, shape, a, b, c, d, row["unit"], row["source"], band, detailed_band, name,
                  row.get("low_text", template.low_text if template else ""),
                  row.get("high_text", template.high_text if template else ""), row["notes"], row["factor"])


def _unresolved(row: dict) -> str | None:
    if any(value.upper() == "TBD" for value in row.values()):
        return "TBD value"
    if not row["source"] or "PROVISIONAL" in row["source"].upper() or "PROVISIONAL" in row.get("low_text", "").upper():
        return "missing or PROVISIONAL source"
    needed = ("a", "b", "c", "d") if row["shape"] == "trapezoid" else ("a", "b")
    if any(not row[key] for key in (*needed, "unit")):
        return "incomplete threshold or unit"
    return None


def load_config(*, soil_path: Path, crops_path: Path, requirements_path: Path, params_path: Path, manifest: dict) -> ScoringConfig:
    paths = {"soil_thresholds": Path(soil_path), "crops": Path(crops_path),
             "crop_requirements": Path(requirements_path), "params": Path(params_path)}
    params = validate_params(yaml.safe_load(paths["params"].read_text()))
    records = manifest.get("bands", [])
    units = {}
    for record in records:
        if not isinstance(record, dict) or record.get("name") in units or not record.get("name"):
            raise ConfigurationError("manifest has missing or duplicate band names")
        units[record["name"]] = canonical_unit(record.get("unit", ""))
    baseline = {"ethiopia_mask", "landcover", "es_confidence", "kiremt_mean_mm", "belg_mean_mm", "kiremt_cv", "belg_cv"}
    if not baseline <= set(units):
        raise ConfigurationError("manifest lacks required mask/land-cover/confidence/rainfall bands")
    soil_rows = _rows(paths["soil_thresholds"], SOIL_COLUMNS)
    factors = []
    skipped = []
    seen = set()
    policies = set(params["active_soil_factors"]) | set(params["skipped_soil_factors"])
    for row in soil_rows:
        name = row["factor"]
        if name in seen or name not in SOIL_BINDINGS:
            raise ConfigurationError(f"duplicate or unknown soil factor: {name}")
        seen.add(name)
        if name not in policies:
            raise ConfigurationError(f"{name}: explicit activation/skip policy is missing")
        reason = params["skipped_soil_factors"].get(name) or _unresolved(row)
        if not reason and ("unconfirmed" in row["notes"].lower() or "recommend not scoring" in row["notes"].lower() or "skip in" in row["notes"].lower()):
            if name not in params["confirmed_units"]:
                reason = "legacy unit/activation note lacks explicit owner confirmation"
        if reason:
            skipped.append({"kind": "soil_factor", "factor": name, "reason": reason})
            continue
        expected = SOIL_BINDINGS[name]
        if bool(row["band"]) == bool(row["detailed_band"]) or (row["band"] or row["detailed_band"]) != expected:
            raise ConfigurationError(f"{name}: must bind only to {expected}")
        if expected not in units:
            raise ConfigurationError(f"{name}: manifest lacks {expected}")
        confirmed = params["confirmed_units"].get(name, row["unit"])
        if canonical_unit(row["unit"]) != units[expected] or canonical_unit(confirmed) != units[expected]:
            skipped.append({"kind": "soil_factor", "factor": name, "reason": "unit incompatible with stack/owner confirmation"})
            continue
        factors.append(_factor(row, name, band=row["band"], detailed_band=row["detailed_band"]))
    if policies - seen:
        raise ConfigurationError("activation policy refers to absent threshold rows")
    if not any(f.band for f in factors) or not any(f.factor == "slope" for f in factors):
        raise ConfigurationError("at least one basic factor and an active slope rule are required")
    generic = {f.factor: f for f in factors}
    all_crops = []
    crop_ids = set()
    for index, row in enumerate(_rows(paths["crops"], {"crop", "display_name", "seasons"}), start=1):
        name = row["crop"]
        seasons = tuple(row["seasons"].split(";"))
        if not re.fullmatch(r"[a-z][a-z0-9_]*", name) or name in crop_ids or index >= 255 or not row["display_name"]:
            raise ConfigurationError("invalid, duplicate or excessive crop IDs")
        if len(seasons) != len(set(seasons)) or not set(seasons) <= set(SEASONS):
            raise ConfigurationError(f"{name}: crop seasons must be unique kiremt/belg values")
        crop_ids.add(name)
        all_crops.append(Crop(name, row["display_name"], seasons, index))
    requirements = {crop.crop: {} for crop in all_crops}
    seen_pairs = set()
    for row in _rows(paths["crop_requirements"], CROP_COLUMNS):
        crop = row["crop"]
        name = params["crop_factor_aliases"].get(row["factor"], row["factor"])
        pair = (crop, name)
        if crop not in crop_ids or name not in set(GATE_ORDER) | (set(SOIL_BINDINGS) - {"slope", "coarse_fragments"}) | {"slope_deg"} or pair in seen_pairs:
            raise ConfigurationError(f"unknown or duplicate crop/factor: {pair}")
        seen_pairs.add(pair)
        reason = _unresolved(row)
        if reason:
            skipped.append({"kind": "crop_factor", "crop": crop, "factor": name, "reason": reason})
            continue
        template = generic.get("slope" if name == "slope_deg" else name)
        if name not in GATE_ORDER and template is None:
            skipped.append({"kind": "crop_factor", "crop": crop, "factor": name, "reason": "generic factor is inactive"})
            continue
        unit = canonical_unit(row["unit"])
        if name == "rain_reliability":
            percent = unit == "%" and row["shape"] == "ramp_up" and params["rainfall_reliability"]["method"] == "normal_fraction_of_mean"
            cv = unit == "fraction" and row["shape"] == "ramp_down"
            if not (percent or cv):
                raise ConfigurationError(f"{crop}: reliability shape/unit incompatible with configured CV-derived representation")
            band = "<season>_cv"
        elif name == "rain_season_mm":
            band = "<season>_mean_mm"
        elif name == "temperature_season_c":
            band = "tmean_<season>_c"
        elif name in {"elevation_m", "clay_pct", "ph"}:
            band = {"elevation_m": "elevation_m", "clay_pct": "sg_clay", "ph": "sg_phh2o"}[name]
        else:
            band = template.input_band
        if name != "rain_reliability":
            check_bands = [band.replace("<season>", season) for season in SEASONS] if "<season>" in band else [band]
            if any(check not in units or units[check] != unit for check in check_bands):
                raise ConfigurationError(f"{crop}/{name}: unit or input band mismatch")
        factor = _factor(row, name, band="" if template and template.detailed_band else band,
                         detailed_band=band if template and template.detailed_band else "", template=template)
        requirements[crop][name] = factor
    computed = []
    for crop in all_crops:
        if not {"elevation_m", "rain_season_mm"} <= set(requirements[crop.crop]):
            skipped.append({"kind": "crop", "crop": crop.crop, "reason": "missing active elevation_m/rain_season_mm requirements"})
            continue
        computed.append(crop)
    inputs = {key: {"file": str(path), "sha256": sha256(path)} for key, path in paths.items()}
    return ScoringConfig(tuple(factors), tuple(computed), requirements, params, units, tuple(skipped), inputs)
