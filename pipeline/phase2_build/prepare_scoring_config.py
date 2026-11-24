#!/usr/bin/env python3
"""Validate and materialize an explicitly approved YAML; never alter the CSVs.

The input YAML contains owner-supplied scientific values. This script contains
none. It logs the exact input hashes and active/skipped rules under interim.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scoring.config import load_config, sha256, validate_params  # noqa: E402


def materialize(source: Path, output: Path, *, replace: bool = False) -> None:
    source, output = source.resolve(), output.resolve()
    if output.is_relative_to(ROOT / "data/inputs") or source == output:
        raise ValueError("source inputs cannot be overwritten or materialized in place")
    validate_params(yaml.safe_load(source.read_text()))
    if output.exists() and output.read_bytes() != source.read_bytes() and not replace:
        raise ValueError("existing params differ; review them before explicitly using --replace")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    shutil.copyfile(source, temporary)
    if sha256(temporary) != sha256(source):
        raise ValueError("configuration copy checksum mismatch")
    temporary.replace(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="YAML containing explicitly owner-approved values")
    parser.add_argument("--output", type=Path, default=ROOT / "data/scoring/params.yaml")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/processed/band_manifest.json")
    parser.add_argument("--soil-thresholds", type=Path, default=ROOT / "data/scoring/soil_thresholds.csv")
    parser.add_argument("--crops", type=Path, default=ROOT / "data/crops/crops.csv")
    parser.add_argument("--crop-requirements", type=Path, required=True, help="Explicit path, including the approved source filename")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "data/interim/score_build")
    parser.add_argument("--replace", action="store_true", help="Explicit permission to replace a differing parameter file")
    args = parser.parse_args()
    # Validate the complete source configuration before modifying the target.
    manifest = json.loads(args.manifest.read_text())
    config = load_config(soil_path=args.soil_thresholds, crops_path=args.crops,
                         requirements_path=args.crop_requirements, params_path=args.source, manifest=manifest)
    materialize(args.source, args.output, replace=args.replace)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(__file__).resolve(), args.work_dir / Path(__file__).name)
    report = {"status": "passed", "owner_approval": config.params["approval"],
              "materialized_params": {"file": str(args.output), "sha256": sha256(args.output)},
              "source_inputs": config.inputs, "active_soil_factors": [f.factor for f in config.soil_factors],
              "computed_crops": [crop.crop for crop in config.crops], "skipped": list(config.skipped),
              "crop_factor_aliases": config.params["crop_factor_aliases"],
              "csv_files_unchanged": all(sha256(Path(v["file"])) == v["sha256"] for k, v in config.inputs.items() if k != "params")}
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    (args.work_dir / "configuration_audit.json").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
