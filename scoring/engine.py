"""One evaluator for offline raster windows and future one-cell explanations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config.grid import NODATA
from scoring import crops, soil
from scoring.config import ADDITIONAL_BASIC, GATE_ORDER, Crop, Factor, ScoringConfig, canonical_unit
from scoring.functions import rainfall_reliability_normal


@dataclass
class Evaluation:
    bands: dict[str, np.ndarray]
    main_season: np.ndarray
    soil_factors: dict
    crop_details: dict


def decode(values) -> np.ndarray:
    """Convert physical nodata to computational NaN without rescaling anything."""
    array = np.asarray(values, dtype=np.float32)
    return np.where(np.isfinite(array) & (array != NODATA), array, np.nan).astype(np.float32)


class ScoringEngine:
    def __init__(self, config: ScoringConfig, calibration: dict):
        self.config = config
        self.calibration = calibration
        if not np.isfinite(calibration["cv_low"]) or not np.isfinite(calibration["cv_high"]) or calibration["cv_low"] >= calibration["cv_high"]:
            raise ValueError("engine requires finite ordered data-derived CV cut-offs")
        self.generic = {factor.factor: factor for factor in config.soil_factors}
        self.soil_lookup = {factor.factor: index for index, factor in enumerate(config.soil_factors, start=1)}
        self.crop_lookups = {}
        for crop in config.crops:
            gates, basic, extra = self.crop_rules(crop)
            names = [*gates, *basic, *extra]
            if len(names) >= 255 or len(names) != len(set(names)):
                raise ValueError("limiting-factor IDs must be unique and fit 1..254")
            self.crop_lookups[crop.crop] = {name: index for index, name in enumerate(names, start=1)}

    def crop_rules(self, crop: Crop) -> tuple[dict[str, Factor], dict[str, Factor], dict[str, Factor]]:
        rows = self.config.requirements[crop.crop]
        active = dict(rows)
        explicitly_skipped = {item["factor"] for item in self.config.skipped if item["kind"] == "crop_factor" and item["crop"] == crop.crop}
        if "ph" not in active and "ph" in self.generic and "ph" not in explicitly_skipped:
            active["ph"] = self.generic["ph"]
        if "rain_reliability" not in active and "rain_reliability" not in explicitly_skipped:
            active["rain_reliability"] = Factor("rain_reliability", "ramp_down", self.calibration["cv_low"], self.calibration["cv_high"], None, None,
                                                "fraction", "Selected main-season CV percentiles inside Ethiopia", band="<season>_cv",
                                                display_name="Rainfall reliability", high_text="high variability")
        gates = {name: active[name] for name in GATE_ORDER if name in active}
        basic = {}
        for name in ADDITIONAL_BASIC:
            fallback = self.generic.get("slope" if name == "slope_deg" else name)
            if name in rows:
                basic[name] = rows[name]
            elif fallback is not None and name not in explicitly_skipped:
                basic[name] = fallback
        extra = {factor.factor: rows.get(factor.factor, factor) for factor in self.config.es_factors if factor.factor not in explicitly_skipped}
        if "slope_deg" not in basic:
            raise ValueError(f"{crop.crop}: an active slope rule is required for the hard-zero contract")
        return gates, basic, extra

    @property
    def required_bands(self) -> set[str]:
        names = {"ethiopia_mask", "landcover", "es_confidence", "kiremt_mean_mm", "belg_mean_mm", "kiremt_cv", "belg_cv"}
        names.update(factor.input_band for factor in self.config.soil_factors)
        for crop in self.config.crops:
            groups = self.crop_rules(crop)
            for group in groups:
                for factor in group.values():
                    if "<season>" in factor.input_band:
                        names.update(factor.input_band.replace("<season>", season) for season in crop.seasons)
                    else:
                        names.add(factor.input_band)
        return names

    def _inputs(self, bands: dict) -> tuple[dict[str, np.ndarray], np.ndarray]:
        missing = self.required_bands - set(bands)
        if missing:
            raise ValueError(f"missing named scoring bands: {sorted(missing)}")
        values = {name: decode(bands[name]) for name in self.required_bands}
        shape = values["ethiopia_mask"].shape
        if any(value.shape != shape for value in values.values()):
            raise ValueError("all input bands must have the same shape")
        if not np.all(np.isin(values["ethiopia_mask"], (0, 1))):
            raise ValueError("ethiopia_mask must contain only 0/1")
        inside = values["ethiopia_mask"] == 1
        for name in values.keys() - {"ethiopia_mask"}:
            values[name] = np.where(inside, values[name], np.nan).astype(np.float32)
        return values, inside

    def evaluate(self, bands: dict, *, breakdown: bool = False) -> Evaluation:
        values, inside = self._inputs(bands)
        weights = self.config.params["soil_quality"]
        threshold = self.config.params["blend"]["detailed_label_min_weight"]
        basic_scores = {factor.factor: factor.score(values[factor.input_band]) for factor in self.config.basic_factors}
        es_scores = {factor.factor: factor.score(values[factor.input_band]) for factor in self.config.es_factors}
        quality_basic = soil.soil_quality(list(basic_scores.values()), **weights)
        quality_detailed = soil.soil_quality(list(basic_scores.values()), list(es_scores.values()), **weights) if es_scores else np.full(inside.shape, np.nan, dtype=np.float32)
        weight = soil.data_weight(values["es_confidence"], list(es_scores.values()))
        quality = soil.blend(quality_basic, quality_detailed, weight)
        dominant_detailed = (weight >= threshold) & np.isfinite(quality_detailed)
        basic_limit = soil.limiting_factor(basic_scores, self.soil_lookup, valid=np.isfinite(quality_basic))
        detailed_limit = soil.limiting_factor({**basic_scores, **es_scores}, self.soil_lookup, valid=np.isfinite(quality_detailed))
        limit = np.where(dominant_detailed, detailed_limit, basic_limit).astype(np.float32)
        main_season, cv = soil.main_rainy_season(values["kiremt_mean_mm"], values["belg_mean_mm"], values["kiremt_cv"], values["belg_cv"])
        risk, overall = soil.drought_score(quality, cv, cv_low=self.calibration["cv_low"], cv_high=self.calibration["cv_high"], weight=self.config.params["drought"]["weight"])
        outputs = {"soil_quality_basic": quality_basic, "soil_quality_detailed": quality_detailed, "soil_quality": quality,
                   "soil_limiting_factor": limit, "data_weight": weight, "drought_risk": risk, "soil_score_overall": overall}
        soil_factors = {}
        if breakdown:
            soil_factors = {factor.factor: {"rule": factor, "value": values[factor.input_band],
                                            "score": (basic_scores | es_scores)[factor.factor],
                                            "group": "basic" if factor.band else "detailed"} for factor in self.config.soil_factors}
        details = {}
        viable = {}
        viability_params = self.config.params["viability"]
        landcover_valid = np.isfinite(values["landcover"])
        permitted = inside & landcover_valid & ~np.isin(values["landcover"], (50, 80))
        for crop in self.config.crops:
            gates, additional_basic, additional_es = self.crop_rules(crop)
            rain_rule = gates["rain_season_mm"]
            means = [values[f"{season}_mean_mm"] for season in crop.seasons]
            chosen = crops.choose_season([rain_rule.score(mean) for mean in means], means)
            gate_values = {}
            for name, rule in gates.items():
                if "<season>" in rule.input_band:
                    value = crops.selected_values([values[rule.input_band.replace("<season>", season)] for season in crop.seasons], chosen)
                else:
                    value = values[rule.input_band]
                if name == "rain_reliability" and canonical_unit(rule.unit) == "%":
                    value = rainfall_reliability_normal(value, threshold_ratio=self.config.params["rainfall_reliability"]["threshold_ratio"])
                gate_values[name] = value
            gate_scores = {name: rule.score(gate_values[name]) for name, rule in gates.items()}
            additional_values = {name: values[rule.input_band] for name, rule in (additional_basic | additional_es).items()}
            b_scores = {name: rule.score(additional_values[name]) for name, rule in additional_basic.items()}
            e_scores = {name: rule.score(additional_values[name]) for name, rule in additional_es.items()}
            b_quality = soil.soil_quality(list(b_scores.values()), **weights)
            d_quality = soil.soil_quality(list(b_scores.values()), list(e_scores.values()), **weights) if e_scores else np.full(inside.shape, np.nan, dtype=np.float32)
            additional = soil.blend(b_quality, d_quality, weight)
            crop_score, gate, modifier = crops.viability(list(gate_scores.values()), additional,
                                                       slope_score=b_scores["slope_deg"], floor=viability_params["fertility_floor"])
            valid = permitted & np.isfinite(crop_score)
            crop_score = np.where(valid, crop_score, np.nan).astype(np.float32)
            lookup = self.crop_lookups[crop.crop]
            gate_limit = soil.limiting_factor(gate_scores, lookup, valid=valid)
            b_limit = soil.limiting_factor(b_scores, lookup, valid=valid)
            d_limit = soil.limiting_factor(b_scores | e_scores, lookup, valid=valid & np.isfinite(d_quality))
            crop_detailed = (weight >= threshold) & np.isfinite(d_quality)
            additional_limit = np.where(crop_detailed, d_limit, b_limit)
            crop_limit = np.where(gate <= modifier, gate_limit, additional_limit)
            crop_limit = np.where(b_scores["slope_deg"] == 0, lookup["slope_deg"], crop_limit)
            crop_limit = np.where(valid, crop_limit, np.nan).astype(np.float32)
            viable[crop.crop] = crop_score
            outputs[f"viability_{crop.crop}"] = crop_score
            outputs[f"viability_{crop.crop}_limiting"] = crop_limit
            if breakdown:
                rules = gates | additional_basic | additional_es
                all_values = gate_values | additional_values
                all_scores = gate_scores | b_scores | e_scores
                details[crop.crop] = {"chosen_season": chosen, "seasons": crop.seasons, "gate": gate,
                                      "modifier": modifier, "additional_basic": b_quality, "additional_detailed": d_quality,
                                      "additional_score": additional, "dominant_detailed": crop_detailed,
                                      "factors": {name: {"rule": rule, "value": all_values[name], "score": all_scores[name],
                                                          "group": "gate" if name in gates else "additional"}
                                                  for name, rule in rules.items()}}
        if viable:
            outputs.update(crops.best_crop(viable, {crop.crop: crop.index for crop in self.config.crops}, minimum_score=viability_params["best_crop_min_score"]))
        else:
            outputs.update({"best_crop": np.full(inside.shape, 255, dtype=np.uint8),
                            "best_crop_score": np.full(inside.shape, np.nan, dtype=np.float32),
                            "best_crop_margin": np.full(inside.shape, np.nan, dtype=np.float32)})
        for name in outputs.keys() - {"best_crop"}:
            outputs[name] = np.where(inside, outputs[name], np.nan).astype(np.float32)
        return Evaluation(outputs, main_season, soil_factors, details)

    def score_cell(self, values: dict) -> dict:
        """JSON-ready one-cell breakdown using precisely the window evaluator."""
        if any(np.asarray(value).shape for value in values.values()):
            raise ValueError("score_cell expects scalar band values")
        evaluated = self.evaluate(values, breakdown=True)
        scores = {name: (None if name == "best_crop" and value == 255 else _json_number(value)) for name, value in evaluated.bands.items()}
        weight = scores["data_weight"]
        cutoff = self.config.params["blend"]["detailed_label_min_weight"]
        tier = None if weight is None else "detailed" if weight >= cutoff else "partial" if weight > 0 else "basic"
        soil_limit = scores["soil_limiting_factor"]
        result = {"scores": scores, "soil": {"tier": tier, "data_weight": weight,
                  "limiting_factor": _lookup_name(self.soil_lookup, soil_limit),
                  "factors": [_factor_cell(name, record, self.soil_lookup, soil_limit) for name, record in evaluated.soil_factors.items()]},
                  "main_season": {0: None, 1: "kiremt", 2: "belg"}[int(evaluated.main_season)], "crops": {}}
        for name, detail in evaluated.crop_details.items():
            index = int(detail["chosen_season"])
            limit = scores[f"viability_{name}_limiting"]
            result["crops"][name] = {"score": scores[f"viability_{name}"], "tier": tier,
                                    "chosen_season": detail["seasons"][index] if index >= 0 else None,
                                    "gate": _json_number(detail["gate"]), "modifier": _json_number(detail["modifier"]),
                                    "additional_basic": _json_number(detail["additional_basic"]),
                                    "additional_detailed": _json_number(detail["additional_detailed"]),
                                    "additional_score": _json_number(detail["additional_score"]),
                                    "limiting_factor": _lookup_name(self.crop_lookups[name], limit),
                                    "factors": [_factor_cell(factor, record, self.crop_lookups[name], limit) for factor, record in detail["factors"].items()]}
        return result


def _json_number(value) -> float | None:
    return float(value) if np.isfinite(value) else None


def _lookup_name(lookup: dict[str, int], factor_id) -> str | None:
    return next((name for name, index in lookup.items() if index == factor_id), None)


def _factor_cell(name: str, record: dict, lookup: dict, limit) -> dict:
    rule = record["rule"]
    value = _json_number(record["value"])
    direction = None
    if value is not None:
        if rule.shape == "ramp_up":
            direction = "low" if value < rule.b else "optimal"
        elif rule.shape == "ramp_down":
            direction = "high" if value > rule.a else "optimal"
        else:
            direction = "low" if value < rule.b else "high" if value > rule.c else "optimal"
    text = rule.low_text if direction == "low" else rule.high_text if direction == "high" else None
    return {"id": name, "name": rule.display_name, "value": value, "score": _json_number(record["score"]),
            "unit": rule.unit, "source": rule.source, "band": rule.input_band, "group": record["group"],
            "shape": rule.shape, "zones": {key: getattr(rule, key) for key in ("a", "b", "c", "d")},
            "ideal_min": rule.b if rule.shape != "ramp_down" else None,
            "ideal_max": rule.c if rule.shape == "trapezoid" else rule.a if rule.shape == "ramp_down" else None,
            "limiting": lookup[name] == limit, "direction": direction, "direction_text": text,
            "resolution_label": "Nominal 1 km cell; EthioSIS values are estimates, not exact measured samples"}
