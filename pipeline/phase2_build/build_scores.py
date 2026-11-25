#!/usr/bin/env python3
"""Build a reviewed-input stack-v2 COG candidate, never publish it automatically.

Run with the approved project .venv. All raster I/O uses names and bounded
windows; no native inputs, IDW interpolation or climate layers are regenerated.
The candidate, byte-exact v1 manifest archive, script snapshot and evidence stay
under data/interim/score_build. --validate-only performs read-only acceptance.
Promotion to data/processed requires a separate owner decision on these results.
"""

from __future__ import annotations

import argparse
import copy
import platform
import sys
import zlib
from dataclasses import replace
from pathlib import Path

import numpy as np
import scipy
from osgeo import gdal

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config import stack, stack_v2  # noqa: E402
from pipeline.phase2_build import diagnostics, provenance, raster, technical  # noqa: E402
from scoring.config import load_config  # noqa: E402
from scoring.engine import ScoringEngine  # noqa: E402

gdal.UseExceptions()

# Engineering I/O budget only: retain this grid's multiband decoded COG tiles
# during full-band acceptance instead of repeatedly evicting/decompressing them.
GDAL_CACHE_BYTES = 1024 * 1024 * 1024


def code_files(root: Path) -> tuple[Path, ...]:
    return tuple(sorted({root / "pipeline/grid.yaml", *(root / "config").glob("*.py"),
                         *(root / "scoring").glob("*.py"), *(root / "pipeline/phase2_build").glob("*.py")}))


def verify_fingerprints(recorded: dict, root: Path, description: str) -> None:
    if not isinstance(recorded, dict) or not recorded:
        raise ValueError(f"{description}: missing input fingerprints")
    observed = provenance.fingerprints((root / name for name in recorded), root)
    if observed != recorded:
        changed = sorted(name for name in set(recorded) | set(observed) if recorded.get(name) != observed.get(name))
        raise ValueError(f"{description}: checksum/size mismatch in {changed}")


def accepted_v1(source_path: Path, manifest_path: Path, validation_path: Path, root: Path) -> tuple[dict, dict]:
    manifest = provenance.load_json(manifest_path)
    stack.validate_manifest(manifest)
    if provenance.sha256(source_path) != manifest.get("raster_sha256") or source_path.stat().st_size != manifest.get("raster_size_bytes"):
        raise ValueError("stack v1 no longer matches its accepted raster checksum/size")
    validation = provenance.load_json(validation_path)
    if (validation.get("status") != "passed" or validation.get("stack_version") != stack.STACK_VERSION
            or validation.get("core_bands") != len(stack.CORE_BAND_NAMES)
            or validation.get("raster_sha256") != manifest["raster_sha256"]
            or validation.get("all_original_inputs_unchanged") is not True):
        raise ValueError("matching passed Phase 1 raster acceptance evidence is required")
    verify_fingerprints(manifest.get("build_inputs"), root, "immutable Phase 1 provenance")
    lookup = manifest.get("lookup", {})
    if not lookup.get("file") or provenance.sha256(root / lookup["file"]) != lookup.get("sha256"):
        raise ValueError("the accepted woreda lookup checksum changed")
    return manifest, validation


def build_manifest(v1_manifest: dict, engine: ScoringEngine, inputs: dict, *, source_path: Path, archive_path: Path,
                   validation_path: Path, raster_path: Path, candidate_path: Path, work_dir: Path, root: Path,
                   score_statistics: dict | None = None) -> dict:
    parent = {"stack_version": stack.STACK_VERSION, "raster_file": provenance.label(source_path, root),
              "raster_sha256": v1_manifest["raster_sha256"], "manifest_file": provenance.label(archive_path, root),
              "manifest_sha256": provenance.sha256(archive_path),
              "band_records_sha256": stack_v2.metadata_sha256(v1_manifest["bands"])}
    records = copy.deepcopy(v1_manifest["bands"]) + stack_v2.score_records(engine, parent["raster_file"])
    if score_statistics is not None:
        for record in records[len(stack.CORE_BAND_NAMES):]:
            record["statistics"] = copy.deepcopy(score_statistics[record["name"]])
    manifest = {"stack_version": stack_v2.STACK_VERSION, "status": "candidate_only", "grid": copy.deepcopy(v1_manifest["grid"]),
                "raster_file": provenance.label(candidate_path, root), "raster_sha256": provenance.sha256(raster_path),
                "raster_size_bytes": raster_path.stat().st_size, "overview_factors": list(v1_manifest["overview_factors"]),
                "parent_stack": parent, "bands": records, "scoring": stack_v2.scoring_metadata(engine),
                "storage": copy.deepcopy(v1_manifest["storage"]), "lookup": copy.deepcopy(v1_manifest["lookup"]),
                "build_script": provenance.BUILD_SCRIPT, "build_inputs": inputs,
                "script_snapshot": provenance.label(work_dir / "script_snapshot", root),
                "parent_acceptance": {"file": provenance.label(validation_path, root), "sha256": provenance.sha256(validation_path)},
                "runtime": {"python": platform.python_version(), "gdal": gdal.VersionInfo("RELEASE_NAME"),
                            "numpy": np.__version__, "scipy": scipy.__version__, "zlib": zlib.ZLIB_RUNTIME_VERSION},
                "build_options": {"window_size": stack.BLOCK_SIZE, "threads": 1,
                                  "gdal_cache_bytes": GDAL_CACHE_BYTES,
                                  "diagnostic_read_source": "bitwise_verified_assembly",
                                  "verification_window_size": technical.VERIFY_WINDOW_SIZE,
                                  "random_inside_cells": technical.RANDOM_INSIDE_CELLS, "random_seed": technical.RANDOM_SEED},
                "assumptions": [
                    "Every grid value comes from pipeline/grid.yaml; accepted v1 is read only, with no resampling, crop, interpolation or native regeneration.",
                    "All 38 v1 band records, original full-resolution float32 bits and per-band provenance tags are retained unchanged.",
                    "Thresholds, units, sources, aliases, tuning and original crop indices are frozen from explicitly approved configuration files; skipped rows remain disclosed.",
                    "Rainfall reliability is the approved Normal fraction-of-seasonal-mean approximation in percent, NOT a probability of satisfying crop water demand; raw seasonal CV is unchanged.",
                    f"Phosphorus is AfricanSoilGrids mg/100kg as supplied, sg_n is the sole nitrogen input in mg/kg, and {len(engine.config.es_factors)} active EthioSIS scored nutrients retain mg/kg values; silicon is informational only.",
                    "Crop CSV temperature explicitly maps to temperature_season_c; existing tmean_<season>_c inputs are neither renamed nor re-derived.",
                    "Chickpea uses the owner-approved Kiremt recharge/residual-moisture proxy with its supplied clay gate, not a simulated soil-water balance.",
                    "Basic/gate missing inputs propagate even beside a hard zero. Missing detailed nutrients are ignored only in their tier and reduce blend completeness; there is no imputation.",
                    "Water, built-up and unknown land-cover cells have no crop score; current cropland is not a viability filter.",
                    "Logical categorical missing 255 is encoded as physical -9999. Continuous scores are 0-1; all JSON missing values are null.",
                    "Every average/mode overview is generated directly from full-resolution values; full-resolution country mask/nodata is authoritative at mixed overview edges.",
                    "Statistics and diagnostics read the band-interleaved assembly only after every original/appended full-resolution pixel is verified bitwise against the candidate COG; reported hashes identify that COG.",
                    "Supplied EthioSIS estimates and their 20 km coverage are not regenerated, expanded or shrunk; radius sensitivity is skipped at the owner's request.",
                    "All values describe nominal 1 km advisory estimates, not exact native 90 m samples or agronomic guarantees. Named objects, access bands and Phase 3 remain deferred.",
                    "This is an interim candidate. Numerical acceptance and scientific diagnostic flags do not authorize promotion without owner review of this exact raster/manifest.",
                ], "deferred_optional_bands": copy.deepcopy(v1_manifest.get("deferred_optional_bands", [])),
                "scope": "Phase 2 offline scoring candidate only; no Django/API/UI, PostGIS, PMTiles or deferred access work"}
    stack_v2.validate_manifest(manifest, v1_manifest=v1_manifest)
    return manifest


def reproducibility_check(source, v1_manifest: dict, manifest: dict, inputs: dict, *, source_path: Path, archive_path: Path,
                          validation_path: Path, candidate_path: Path, raw_path: Path, work_dir: Path, root: Path) -> dict:
    engine = stack_v2.frozen_engine(manifest)
    fresh = raster.calibrate_stack(source, stack.band_indices_by_name(source), engine.config.params,
                                   window_size=technical.VERIFY_WINDOW_SIZE)
    if fresh != engine.calibration:
        raise ValueError("reproducibility calibration changed")
    repeat_dir = work_dir / "reproducibility"
    repeat_dir.mkdir(parents=True, exist_ok=True)
    repeated_raw, repeated_cog = repeat_dir / raw_path.name, repeat_dir / candidate_path.name
    raster.write_stack(source, manifest["bands"], engine, repeated_raw, repeated_cog,
                       parent_sha256=manifest["parent_stack"]["raster_sha256"], overview_factors=manifest["overview_factors"])
    repeated = gdal.OpenEx(str(repeated_cog), gdal.OF_RASTER | gdal.OF_READONLY)
    indices = stack.band_indices_by_name(repeated)
    inside = diagnostics.collect_band(repeated, indices, "ethiopia_mask") == 1
    stats = {row["name"]: diagnostics.band_statistics(diagnostics.collect_band(repeated, indices, row["name"]), inside, row)
             for row in manifest["bands"][len(stack.CORE_BAND_NAMES):]}
    second_manifest = build_manifest(v1_manifest, engine, inputs, source_path=source_path, archive_path=archive_path,
        validation_path=validation_path, raster_path=repeated_cog, candidate_path=candidate_path, work_dir=work_dir, root=root,
        score_statistics=stats)
    repeated = None
    first_bytes, second_bytes = provenance.json_bytes(manifest), provenance.json_bytes(second_manifest)
    if provenance.sha256(repeated_cog) != manifest["raster_sha256"] or first_bytes != second_bytes:
        raise ValueError("rerun COG/manifest bytes differ; nondeterminism is not accepted silently")
    first_raw, second_raw = provenance.sha256(raw_path), provenance.sha256(repeated_raw)
    return {"status": "byte_identical", "method": "Fresh window calibration and complete independent assembly, COG translation and statistics; identical logical candidate paths and single-threaded compression",
            "raster_sha256": manifest["raster_sha256"], "repeated_raster_sha256": provenance.sha256(repeated_cog),
            "manifest_bytes_equal": True, "raw_assembly_bytes_equal": first_raw == second_raw,
            "first_raw_sha256": first_raw, "repeated_raw_sha256": second_raw,
            "scope": "Proven for the recorded inputs/code/runtime; cross-version GDAL/compression reproducibility is not claimed",
            "repeated_raster_file": provenance.label(repeated_cog, root)}


def configuration_readiness(engine: ScoringEngine) -> dict:
    intentional = engine.config.params["skipped_soil_factors"]
    unresolved = [row for row in engine.config.skipped
                  if row["kind"] != "soil_factor" or row["factor"] not in intentional or row["reason"] != intentional[row["factor"]]]
    return {"configuration_complete": not unresolved, "unresolved_skipped_rows": unresolved,
            "active_soil_factors": [row.factor for row in engine.config.soil_factors],
            "active_crops": [crop.crop for crop in engine.config.crops],
            "promotion": "blocked pending owner review of actual distributions, sensitivity and diagnostic flags"}


def build(args) -> dict:
    native_before = provenance.input_inventory(ROOT)
    source_manifest_hash = provenance.sha256(args.source_manifest)
    scripts = code_files(ROOT)
    work_dir, managed = provenance.check_workspace(args.work_dir, ROOT, code_files=scripts)
    v1_manifest, references = accepted_v1(args.source_stack, args.source_manifest, args.v1_validation, ROOT)
    config = load_config(soil_path=args.soil_thresholds, crops_path=args.crops, requirements_path=args.crop_requirements,
                         params_path=args.params, manifest=v1_manifest)
    config = replace(config, inputs={key: {**row, "file": provenance.label(ROOT / row["file"], ROOT)} for key, row in config.inputs.items()})
    archive_path = work_dir / "band_manifest_v1.json"
    protected = {args.source_stack, args.source_manifest, args.v1_validation, *(ROOT / name for name in v1_manifest["build_inputs"]),
                 ROOT / v1_manifest["lookup"]["file"], *(ROOT / row["file"] for row in config.inputs.values()), *scripts}
    targets = {(work_dir / name).resolve() for name in provenance.workspace_targets(managed)}
    if any(path.resolve() in targets for path in protected):
        raise ValueError("an output path overlaps a protected source input")
    if archive_path.exists() and archive_path.read_bytes() != args.source_manifest.read_bytes():
        raise ValueError("existing v1 archive differs; do not replace it or its workspace state")
    source = gdal.OpenEx(str(args.source_stack), gdal.OF_RASTER | gdal.OF_READONLY)
    indices = stack.assert_stack_matches_manifest(source, v1_manifest)
    provenance.begin_workspace(work_dir, managed)
    provenance.archive_manifest(args.source_manifest, archive_path)
    provenance.snapshot_code(scripts, work_dir, ROOT)
    paths = {ROOT / name for name in v1_manifest["build_inputs"]} | set(scripts) | {
        args.source_stack, args.v1_validation, ROOT / v1_manifest["lookup"]["file"], archive_path,
        *(ROOT / row["file"] for row in config.inputs.values())}
    inputs = provenance.fingerprints(paths, ROOT)
    for name, expected in v1_manifest["build_inputs"].items():
        if inputs[name] != expected:
            raise ValueError(f"accepted v1 input changed before assembly: {name}")
    source_label = provenance.label(args.source_stack, ROOT)
    if inputs[source_label] != {"sha256": v1_manifest["raster_sha256"], "size_bytes": v1_manifest["raster_size_bytes"]}:
        raise ValueError("accepted v1 raster changed before assembly")
    lookup_label = provenance.label(ROOT / v1_manifest["lookup"]["file"], ROOT)
    if inputs[lookup_label]["sha256"] != v1_manifest["lookup"]["sha256"]:
        raise ValueError("accepted woreda lookup changed before assembly")
    for path in scripts:
        name = provenance.label(path, ROOT)
        if provenance.sha256(work_dir / "script_snapshot" / name) != inputs[name]["sha256"]:
            raise ValueError(f"script changed during snapshotting: {name}")
    for row in config.inputs.values():
        if inputs[row["file"]]["sha256"] != row["sha256"]:
            raise ValueError("configuration input changed while it was being loaded")
    raw_path, candidate_path = work_dir / "common_stack_v2_source.tif", work_dir / "common_stack_v2.tif"
    with (work_dir / "build.log").open("w", encoding="utf-8") as log:
        def note(message):
            print(message, flush=True)
            log.write(message + "\n")
            log.flush()

        note(f"Accepted immutable v1: {len(stack.CORE_BAND_NAMES)} bands; source and audit checksums passed.")
        calibration = raster.calibrate_stack(source, indices, config.params)
        engine = ScoringEngine(config, calibration)
        note(f"Active scoring: {len(config.crops)} crops, {len(config.soil_factors)} soil factors; {len(config.skipped)} disclosed skips.")
        note(f"CV calibration: low={calibration['cv_low']:.9g}, high={calibration['cv_high']:.9g}; {calibration['valid_calibration_cells']:,} valid inside-country cells.")
        records = copy.deepcopy(v1_manifest["bands"]) + stack_v2.score_records(engine, provenance.label(args.source_stack, ROOT))
        note("Assembling candidate and direct-from-full-resolution average/mode COG overviews.")
        raster.write_stack(source, records, engine, raw_path, candidate_path,
                           parent_sha256=v1_manifest["raster_sha256"], overview_factors=v1_manifest["overview_factors"])
        manifest = build_manifest(v1_manifest, engine, inputs, source_path=args.source_stack, archive_path=archive_path,
            validation_path=args.v1_validation, raster_path=candidate_path, candidate_path=candidate_path, work_dir=work_dir, root=ROOT)
        candidate = gdal.OpenEx(str(candidate_path), gdal.OF_RASTER | gdal.OF_READONLY)
        raw = gdal.OpenEx(str(raw_path), gdal.OF_RASTER | gdal.OF_READONLY)
        note("Checking every full-resolution pixel, original v1 tags, COG structure, every overview and scalar reports.")
        acceptance = technical.validate_package(source, candidate, manifest, v1_manifest, raw=raw)
        for row in manifest["bands"][len(stack.CORE_BAND_NAMES):]:
            row["statistics"] = acceptance["appended_band_statistics"][row["name"]]
        stack_v2.validate_manifest(manifest, v1_manifest=v1_manifest)
        provenance.write_json(work_dir / "band_manifest.json", manifest)
        manifest_hash = provenance.sha256(work_dir / "band_manifest.json")
        acceptance.update({"raster_sha256": manifest["raster_sha256"], "manifest_sha256": manifest_hash})
        # Keep useful numerical evidence after an interruption, but never label
        # an incomplete build as accepted or ready for owner review/promotion.
        provenance.write_json(work_dir / "technical_validation.json", {
            **acceptance, "status": "passed_numerical_checks_pending_build_completion",
            "pending": ["scientific diagnostics", "independent reproducibility", "final input-preservation checks"],
        })
        note(f"Technical acceptance passed: {len(manifest['bands'])} bands; {acceptance['scalar_checks']['compared_cells']} scalar checks; no bitwise mismatches.")
        note("Computing diagnostics and approved sensitivity variants from the bitwise-verified band-interleaved assembly.")
        # Both its v1 prefix and score bands have been compared to source/COG
        # above. No source interpolation, scientific formula or values change.
        review = diagnostics.review_diagnostics(raw, raw, manifest, references=references)
        for name, stats in acceptance["appended_band_statistics"].items():
            if review["band_statistics"][name] != stats:
                raise ValueError(f"{name}: independent diagnostic statistics changed")
        candidate = raw = None
        note("Rebuilding independently to require byte-identical candidate raster and manifest.")
        repeated = reproducibility_check(source, v1_manifest, manifest, inputs, source_path=args.source_stack,
            archive_path=archive_path, validation_path=args.v1_validation, candidate_path=candidate_path,
            raw_path=raw_path, work_dir=work_dir, root=ROOT)
        source = None
        verify_fingerprints(inputs, ROOT, "Phase 2 input preservation")
        if provenance.input_inventory(ROOT) != native_before or provenance.sha256(args.source_manifest) != source_manifest_hash:
            raise ValueError("an original input inventory entry or the shared v1 manifest changed during the build")
        if provenance.sha256(work_dir / "band_manifest.json") != manifest_hash:
            raise ValueError("candidate manifest changed during diagnostics/reproducibility")
        report = {"status": "awaiting_owner_review", "stack_version": stack_v2.STACK_VERSION,
                  "raster_file": manifest["raster_file"], "raster_sha256": manifest["raster_sha256"], "manifest_sha256": manifest_hash,
                  "technical": acceptance, "diagnostics": review, "reproducibility": repeated,
                  "calibration": calibration, "skipped": list(config.skipped), "readiness": configuration_readiness(engine),
                  "configuration": stack_v2.scoring_metadata(engine)["configuration"],
                  "preservation": {"input_checksums_unchanged": True, "checksum_verified_files": len(inputs),
                                   "native_input_inventory_unchanged": True, "native_input_files": len(native_before),
                                   "method": "SHA-256/size checks before and after for all v1-provenance and Phase 2 consumed inputs/code; inventory/size/mtime for unused native files, which are never read by scoring",
                                   "v1_raster_unchanged": True, "shared_manifest_unchanged": True,
                                   "source_manifest_file": provenance.label(args.source_manifest, ROOT), "source_manifest_sha256": source_manifest_hash,
                                   "v1_manifest_archive": manifest["parent_stack"]["manifest_file"], "v1_manifest_byte_identical": True},
                  "scope": manifest["scope"], "promotion": "Not performed; data/processed and the shared stack-v1 manifest are unchanged"}
        provenance.write_json(work_dir / "technical_validation.json", acceptance)
        provenance.write_json(work_dir / "score_build_report.json", report)
        provenance.write_bytes(work_dir / "score_build_report.md", diagnostics.markdown_report(report).encode("utf-8"))
        note(f"Reproducibility passed. Review flags: {len(review['warnings'])}. Candidate retained; no promotion performed.")
        note(f"Review report: {provenance.label(work_dir / 'score_build_report.md', ROOT)}")
    provenance.finish_workspace(work_dir, managed)
    return report


def validate_only(work_dir: Path) -> dict:
    work_dir, _ = provenance.check_workspace(work_dir, ROOT, code_files=code_files(ROOT))
    state = provenance.load_json(work_dir / provenance.STATE_FILE)
    if state.get("status") != "awaiting_owner_review":
        raise ValueError("read-only acceptance requires a completed candidate and review report")
    manifest_path = work_dir / "band_manifest.json"
    manifest = provenance.load_json(manifest_path)
    stack_v2.validate_manifest(manifest)
    parent = manifest["parent_stack"]
    archive_path = ROOT / parent["manifest_file"]
    if provenance.sha256(archive_path) != parent["manifest_sha256"]:
        raise ValueError("v1 archive bytes no longer match the frozen parent fingerprint")
    candidate_path, source_path = ROOT / manifest["raster_file"], ROOT / parent["raster_file"]
    if provenance.sha256(candidate_path) != manifest["raster_sha256"] or candidate_path.stat().st_size != manifest["raster_size_bytes"]:
        raise ValueError("candidate raster bytes no longer match its frozen manifest")
    verify_fingerprints(manifest["build_inputs"], ROOT, "candidate provenance")
    v1_manifest, _ = accepted_v1(source_path, archive_path, ROOT / manifest["parent_acceptance"]["file"], ROOT)
    source = gdal.OpenEx(str(source_path), gdal.OF_RASTER | gdal.OF_READONLY)
    candidate = gdal.OpenEx(str(candidate_path), gdal.OF_RASTER | gdal.OF_READONLY)
    raw_path = work_dir / "common_stack_v2_source.tif"
    raw = gdal.OpenEx(str(raw_path), gdal.OF_RASTER | gdal.OF_READONLY) if raw_path.is_file() else None
    acceptance = technical.validate_package(source, candidate, manifest, v1_manifest, raw=raw)
    source = candidate = raw = None
    report = provenance.load_json(work_dir / "score_build_report.json")
    if report.get("raster_sha256") != manifest["raster_sha256"] or report.get("manifest_sha256") != provenance.sha256(manifest_path):
        raise ValueError("review evidence does not identify this exact candidate")
    return {"status": "passed_read_only", "band_count": acceptance["band_count"],
            "scalar_checks": acceptance["scalar_checks"]["compared_cells"], "raster_sha256": manifest["raster_sha256"],
            "manifest_sha256": provenance.sha256(manifest_path), "review_flags": report["diagnostics"]["warnings"],
            "promotion": "not performed"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-stack", type=Path, default=ROOT / "data/processed/common_stack_v1.tif")
    parser.add_argument("--source-manifest", type=Path, default=ROOT / "data/processed/band_manifest.json")
    parser.add_argument("--v1-validation", type=Path, default=ROOT / "data/processed/stack_v1_validation.json")
    parser.add_argument("--soil-thresholds", type=Path, default=ROOT / "data/scoring/soil_thresholds.csv")
    parser.add_argument("--crops", type=Path, default=ROOT / "data/crops/crops.csv")
    parser.add_argument("--crop-requirements", type=Path, default=ROOT / "data/crops/crop_requirments.csv",
                        help="Owner-approved filename on disk; the CSV remains unchanged")
    parser.add_argument("--params", type=Path, default=ROOT / "data/scoring/params.yaml")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "data/interim/score_build")
    parser.add_argument("--validate-only", action="store_true", help="Read-only technical acceptance of the existing candidate")
    args = parser.parse_args()
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, value.resolve())
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetCacheMax(GDAL_CACHE_BYTES)
    try:
        result = validate_only(args.work_dir) if args.validate_only else build(args)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Phase 2 acceptance failed: {exc}\nNo promotion was performed.\n")
    if args.validate_only:
        print(provenance.json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
