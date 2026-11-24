"""Deliberately synthetic ranges and values used only by tests."""

import csv
from pathlib import Path

import yaml

from scoring.config import CROP_COLUMNS, SOIL_COLUMNS, load_config


def params_fixture():
    return {
        "soil_quality": {"mean_weight": 0.5, "min_weight": 0.5},
        "drought": {"weight": 0.2, "cv_low_percentile": 25, "cv_high_percentile": 90, "percentile_method": "linear"},
        "blend": {"detailed_label_min_weight": 0.5},
        "viability": {"fertility_floor": 0.8, "best_crop_min_score": 0.25, "best_crop_tie_margin": 0.05},
        "score_classes": [25, 50, 70, 85],
        "rainfall_reliability": {"method": "normal_fraction_of_mean", "threshold_ratio": 0.75},
        "crop_factor_aliases": {"temperature": "temperature_season_c"},
        "active_soil_factors": ["ph", "organic_carbon", "slope", "zinc"],
        "skipped_soil_factors": {}, "confirmed_units": {"zinc": "mg/kg"},
        "sensitivity": {"mean_min_weights": [[0.3, 0.7]], "fertility_floors": [0.7], "drought_weights": [0.0], "coverage_radius": "skip_owner_requested"},
        "approval": {"approved": True, "authority": "Synthetic test fixture only", "chickpea_rainfall_proxy": "Test fixture"},
    }


def manifest_fixture():
    units = {"ethiopia_mask": "0 or 1", "landcover": "class code", "es_confidence": "fraction",
             "sg_phh2o": "pH", "sg_soc": "g/kg", "slope_mean_deg": "degrees", "es_zn": "mg/kg",
             "elevation_m": "m", "sg_clay": "%"}
    for season in ("kiremt", "belg"):
        units.update({f"{season}_mean_mm": "mm", f"{season}_cv": "fraction", f"tmean_{season}_c": "degrees Celsius"})
    return {"bands": [{"name": name, "unit": unit} for name, unit in units.items()]}


def soil_rows_fixture():
    result = []
    for name, band, detailed, shape, a, b, c, d, unit in (
        ("ph", "sg_phh2o", "", "trapezoid", 3, 5, 7, 9, "pH"),
        ("organic_carbon", "sg_soc", "", "ramp_up", 0, 20, "", "", "g/kg"),
        ("slope", "slope_mean_deg", "", "ramp_down", 5, 30, "", "", "degrees"),
        ("zinc", "", "es_zn", "ramp_up", 0, 2, "", "", "mg/kg"),
    ):
        result.append(dict(factor=name, display_name=name, band=band, detailed_band=detailed, shape=shape,
                           a=a, b=b, c=c, d=d, unit=unit, source="Synthetic test source", notes="",
                           low_text="low", high_text="high"))
    return result


def requirements_fixture(crops=("teff", "wheat")):
    result = []
    for crop in crops:
        for name, shape, a, b, c, d, unit in (
            ("elevation_m", "trapezoid", 0, 1000, 2000, 4000, "m"),
            ("rain_season_mm", "trapezoid", 100, 300, 700, 1000, "mm"),
            ("temperature", "trapezoid", 0, 10, 25, 40, "°C"),
            ("rain_reliability", "ramp_up", 40, 80, "", "", "%"),
            ("ph", "trapezoid", 3, 5, 7, 9, "pH"),
            ("clay_pct", "trapezoid", 0, 10, 50, 80, "%"),
        ):
            result.append(dict(crop=crop, factor=name, shape=shape, a=a, b=b, c=c, d=d,
                               unit=unit, source="Synthetic test source", notes=""))
    return result


def write_csv(path: Path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(columns))
        writer.writeheader()
        writer.writerows(rows)


def write_fixture(root: Path, *, params=None, soil=None, requirements=None):
    paths = {"soil_path": root / "soil.csv", "crops_path": root / "crops.csv",
             "requirements_path": root / "requirements.csv", "params_path": root / "params.yaml"}
    write_csv(paths["soil_path"], SOIL_COLUMNS, soil if soil is not None else soil_rows_fixture())
    write_csv(paths["crops_path"], {"crop", "display_name", "seasons"},
              [{"crop": name, "display_name": name.title(), "seasons": "kiremt;belg"} for name in ("teff", "wheat")])
    write_csv(paths["requirements_path"], CROP_COLUMNS, requirements if requirements is not None else requirements_fixture())
    paths["params_path"].write_text(yaml.safe_dump(params if params is not None else params_fixture()))
    return paths


def config_fixture(root: Path):
    return load_config(**write_fixture(root), manifest=manifest_fixture())


def values_fixture():
    return {"ethiopia_mask": 1, "landcover": 40, "es_confidence": 0.8,
            "sg_phh2o": 6, "sg_soc": 15, "slope_mean_deg": 10, "es_zn": 1,
            "elevation_m": 1500, "sg_clay": 30,
            "kiremt_mean_mm": 600, "belg_mean_mm": 200, "kiremt_cv": 0.2, "belg_cv": 0.5,
            "tmean_kiremt_c": 20, "tmean_belg_c": 15}
