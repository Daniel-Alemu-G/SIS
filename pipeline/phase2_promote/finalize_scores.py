#!/usr/bin/env python3
"""Publish the exact owner-approved stack-v2 candidate without rebuilding rasters.

Approval names both candidate hashes and accepts every disclosed review flag.
The original v1 manifest is archived byte-for-byte before the shared manifest
advances. All raster reads are read-only; copies retain the approved COG bytes.
Generated metadata, logs and publisher snapshots stay under data/interim first.
"""

from __future__ import annotations

import argparse
import copy
import datetime
import hashlib
import shutil
import sys
from pathlib import Path

from osgeo import gdal

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from config import stack, stack_v2  # noqa: E402
from pipeline.phase2_build import build_scores, diagnostics, provenance, technical  # noqa: E402

gdal.UseExceptions()

PUBLISHER = "pipeline/phase2_promote/finalize_scores.py"
STATE_FILE = "publication_state.json"
SCOPE = "Phase 2 scoring stack v2 finalization only"
DELIVERABLES = (
    "band_manifest_v1.json", "common_stack_v2.tif", "band_manifest_v2.json",
    "stack_v2_approval.json", "stack_v2_validation.json", "score_build_report.json",
    "score_build_report.md", "stack_v2_promotion.json", "band_manifest.json",
)
STAGED = ("band_manifest_v2.json", "stack_v2_validation.json", "score_build_report.json",
          "score_build_report.md", "stack_v2_promotion.json", "finalize.log")


def fingerprint(path: Path) -> dict:
    return {"sha256": provenance.sha256(path), "size_bytes": path.stat().st_size}


def bytes_fingerprint(value: bytes) -> dict:
    return {"sha256": hashlib.sha256(value).hexdigest(), "size_bytes": len(value)}


def publisher_files(root: Path) -> tuple[Path, ...]:
    return tuple(sorted((root / "pipeline/phase2_promote").glob("*.py")))


def guard(path: Path, root: Path) -> None:
    """Reject aliases/special files before any narrow, explicitly owned write."""
    if (not path.is_relative_to(root) or path == root or ".." in path.parts
            or not path.resolve().is_relative_to(root.resolve())):
        raise ValueError(f"publication target is outside the project: {path}")
    for parent in (path, *path.parents[:len(path.relative_to(root).parts) - 1]):
        if parent.is_symlink():
            raise ValueError(f"refusing a symlink publication target: {parent}")
    if path.exists() and (not path.is_file() or path.stat().st_nlink != 1):
        raise ValueError(f"publication target is not a regular, singly-linked file: {path}")


def check_approval(approval: dict, manifest: dict, manifest_hash: str, report: dict) -> None:
    required = {"approved", "authority", "date", "statement", "scope",
                "candidate_raster_sha256", "candidate_manifest_sha256", "accepted_review_flags"}
    if not isinstance(approval, dict) or set(approval) != required or approval["approved"] is not True:
        raise ValueError("explicit owner approval of the exact candidate is required")
    for key in ("authority", "date", "statement"):
        if not isinstance(approval[key], str) or not approval[key].strip():
            raise ValueError(f"owner approval requires {key}")
    try:
        datetime.date.fromisoformat(approval["date"])
    except ValueError as exc:
        raise ValueError("owner approval date must be ISO YYYY-MM-DD") from exc
    if approval["scope"] != SCOPE:
        raise ValueError("owner approval must be limited to Phase 2 finalization")
    if (approval["candidate_raster_sha256"] != manifest["raster_sha256"]
            or approval["candidate_manifest_sha256"] != manifest_hash):
        raise ValueError("owner approval does not identify this exact candidate raster/manifest")
    flags = [row["code"] for row in report["diagnostics"]["warnings"]]
    accepted = approval["accepted_review_flags"]
    if not isinstance(accepted, list) or accepted != flags:
        raise ValueError("owner approval must explicitly accept every disclosed review flag in report order")


def candidate_package(candidate_dir: Path, approval_path: Path) -> tuple[dict, dict, dict, dict]:
    candidate_dir, _ = provenance.check_workspace(candidate_dir, ROOT, code_files=build_scores.code_files(ROOT))
    state = provenance.load_json(candidate_dir / provenance.STATE_FILE)
    required_evidence = {"common_stack_v2.tif", "band_manifest.json", "band_manifest_v1.json",
                         "technical_validation.json", "score_build_report.json", "score_build_report.md",
                         "reproducibility/common_stack_v2.tif"}
    if state.get("status") != "awaiting_owner_review" or not required_evidence <= set(state.get("completed_outputs", {})):
        raise ValueError("finalization requires a completed candidate and review evidence")
    manifest_path = candidate_dir / "band_manifest.json"
    manifest, report = provenance.load_json(manifest_path), provenance.load_json(candidate_dir / "score_build_report.json")
    manifest_hash = provenance.sha256(manifest_path)
    stack_v2.validate_manifest(manifest)
    expected_file = provenance.label(candidate_dir / "common_stack_v2.tif", ROOT)
    if manifest.get("status") != "candidate_only" or manifest["raster_file"] != expected_file:
        raise ValueError("only the named completed interim candidate may be promoted")
    if (report.get("status") != "awaiting_owner_review" or report.get("raster_sha256") != manifest["raster_sha256"]
            or report.get("manifest_sha256") != manifest_hash or report.get("technical", {}).get("status") != "passed"
            or report["technical"].get("scalar_checks", {}).get("compared_cells", 0) < 20
            or report.get("reproducibility", {}).get("status") != "byte_identical"
            or report["reproducibility"].get("manifest_bytes_equal") is not True
            or report.get("readiness", {}).get("configuration_complete") is not True
            or report["readiness"].get("unresolved_skipped_rows") != []):
        raise ValueError("completed technical/reproducibility evidence and resolved scoring inputs are required")
    readiness = build_scores.configuration_readiness(stack_v2.frozen_engine(manifest))
    if not readiness["configuration_complete"]:
        raise ValueError("unresolved or provisional rows must not reach publication")
    build_scores.verify_fingerprints(manifest["build_inputs"], ROOT, "approved candidate input preservation")
    if fingerprint(ROOT / manifest["raster_file"]) != {"sha256": manifest["raster_sha256"], "size_bytes": manifest["raster_size_bytes"]}:
        raise ValueError("candidate COG bytes changed after review")
    v1, _ = build_scores.accepted_v1(ROOT / manifest["parent_stack"]["raster_file"],
        ROOT / manifest["parent_stack"]["manifest_file"], ROOT / manifest["parent_acceptance"]["file"], ROOT)
    if provenance.sha256(ROOT / manifest["parent_stack"]["manifest_file"]) != manifest["parent_stack"]["manifest_sha256"]:
        raise ValueError("original v1 manifest archive changed")
    stack_v2.validate_manifest(manifest, v1_manifest=v1)
    approval = provenance.load_json(approval_path)
    check_approval(approval, manifest, manifest_hash, report)
    return manifest, report, v1, approval


def published_manifest(candidate: dict, approval: dict, approval_path: Path, candidate_dir: Path,
                       work_dir: Path, processed: Path, scripts: tuple[Path, ...]) -> dict:
    result = copy.deepcopy(candidate)
    result["status"] = "accepted"
    result["raster_file"] = provenance.label(processed / "common_stack_v2.tif", ROOT)
    result["parent_stack"]["manifest_file"] = provenance.label(processed / "band_manifest_v1.json", ROOT)
    result["scope"] = "Phase 2 accepted offline scoring stack; no Django/API/UI, PostGIS, PMTiles or deferred access work"
    pending = "This is an interim candidate. Numerical acceptance and scientific diagnostic flags do not authorize promotion without owner review of this exact raster/manifest."
    result["assumptions"] = [line for line in result["assumptions"] if line != pending]
    result["assumptions"].append("Owner approved this exact advisory candidate and its disclosed scientific review flags; no threshold, score, raster byte or scoring lookup was altered during publication. This is not independent field/agronomic validation.")
    result["publication"] = {
        "approved": True, "authority": approval["authority"], "date": approval["date"],
        "approval_file": provenance.label(processed / "stack_v2_approval.json", ROOT),
        "approval_sha256": provenance.sha256(approval_path), "accepted_review_flags": approval["accepted_review_flags"],
        "reviewed_candidate": {"raster_file": candidate["raster_file"], "raster_sha256": candidate["raster_sha256"],
            "manifest_file": provenance.label(candidate_dir / "band_manifest.json", ROOT),
            "manifest_sha256": provenance.sha256(candidate_dir / "band_manifest.json"),
            "report_file": provenance.label(candidate_dir / "score_build_report.json", ROOT),
            "report_sha256": provenance.sha256(candidate_dir / "score_build_report.json")},
        "publisher": PUBLISHER, "publisher_inputs": provenance.fingerprints(scripts, ROOT),
        "script_snapshot": provenance.label(work_dir / "script_snapshot", ROOT),
        "workspace": provenance.label(work_dir, ROOT),
        "validation_file": provenance.label(processed / "stack_v2_validation.json", ROOT),
        "receipt_file": provenance.label(processed / "stack_v2_promotion.json", ROOT),
        "manifest_changes": "Publication status/paths, stable v1 archive reference and owner approval only; bands, scoring configuration/calibration, input provenance and raster bytes unchanged",
        "review_meaning": "Owner acceptance with retained advisory limitations and diagnostic flags; no independent agronomic or visual inspection is claimed by the agent",
    }
    stack_v2.validate_manifest(result)
    return result


class PublicationWorkspace:
    """Narrow, checksum-guarded ownership with checkpoints after each write."""

    def __init__(self, work_dir: Path, processed: Path, identity: dict, scripts: tuple[Path, ...]):
        self.work_dir, self.processed = work_dir, processed
        self.state_path = work_dir / STATE_FILE
        snapshots = tuple(work_dir / "script_snapshot" / provenance.label(path, ROOT) for path in scripts)
        self.targets = tuple(processed / name for name in DELIVERABLES) + tuple(work_dir / name for name in STAGED) + snapshots
        checked = (*self.targets, self.state_path, *(path.with_name(path.name + ".tmp") for path in (*self.targets, self.state_path)))
        for path in checked:
            guard(path, ROOT)
        self.state = {"publisher": PUBLISHER, "identity": identity, "status": "finalizing", "completed_outputs": {}}
        if self.state_path.exists():
            previous = provenance.load_json(self.state_path)
            if (not isinstance(previous, dict) or set(previous) != set(self.state)
                    or previous.get("publisher") != PUBLISHER or previous.get("identity") != identity
                    or previous.get("status") not in {"finalizing", "accepted"}
                    or not isinstance(previous.get("completed_outputs"), dict)):
                raise ValueError("publication workspace belongs to an unfamiliar candidate or approval")
            allowed = {provenance.label(path, ROOT) for path in self.targets}
            for name, expected in previous["completed_outputs"].items():
                if name not in allowed or not (ROOT / name).is_file() or fingerprint(ROOT / name) != expected:
                    raise ValueError(f"previous publication output changed outside this step: {name}")
            self.state = previous
        else:
            exceptions = {processed / "band_manifest.json", processed / "band_manifest_v1.json"}
            if any(path.exists() for path in checked if path not in exceptions):
                raise ValueError("existing publication outputs have no known owner; do not overwrite them")

    def checkpoint(self, path: Path | None = None) -> None:
        if path is not None:
            self.state["completed_outputs"][provenance.label(path, ROOT)] = fingerprint(path)
        provenance.write_json(self.state_path, self.state)

    def write(self, path: Path, payload: bytes | Path, *, replace_v1: bool = False) -> None:
        if path not in self.targets or (replace_v1 and path != self.processed / "band_manifest.json"):
            raise ValueError("attempt to write an unowned publication target")
        expected = fingerprint(payload) if isinstance(payload, Path) else bytes_fingerprint(payload)
        if path.exists() and fingerprint(path) == expected:
            self.checkpoint(path)
            return
        if path.exists() and not replace_v1:
            raise ValueError(f"existing publication output differs; do not overwrite it: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        guard(path, ROOT)
        guard(temporary, ROOT)
        if isinstance(payload, Path):
            shutil.copyfile(payload, temporary)
        else:
            temporary.write_bytes(payload)
        if fingerprint(temporary) != expected:
            raise ValueError(f"publication copy is not byte-identical: {path}")
        temporary.replace(path)
        self.checkpoint(path)


def full_validation(manifest: dict, v1: dict, processed: Path) -> dict:
    raster_path = processed / "common_stack_v2.tif"
    expected = {"sha256": manifest["raster_sha256"], "size_bytes": manifest["raster_size_bytes"]}
    if fingerprint(raster_path) != expected:
        raise ValueError("published COG bytes differ from the owner-approved candidate")
    source = gdal.OpenEx(str(ROOT / manifest["parent_stack"]["raster_file"]), gdal.OF_RASTER | gdal.OF_READONLY)
    output = gdal.OpenEx(str(raster_path), gdal.OF_RASTER | gdal.OF_READONLY)
    result = technical.validate_package(source, output, manifest, v1)
    source = output = None
    if fingerprint(raster_path) != expected:
        raise ValueError("published COG bytes changed during read-only acceptance")
    return {**result, "raster_file": manifest["raster_file"], "raster_sha256": manifest["raster_sha256"],
            "manifest_sha256": bytes_fingerprint(provenance.json_bytes(manifest))["sha256"]}


def accepted_report(original: dict, manifest: dict, approval: dict, acceptance: dict) -> dict:
    result = copy.deepcopy(original)
    result.update({"status": "accepted_owner_review", "raster_file": manifest["raster_file"],
                   "manifest_sha256": acceptance["manifest_sha256"], "technical": acceptance,
                   "publication": manifest["publication"], "owner_approval": approval, "scope": manifest["scope"]})
    result["candidate_build_preservation"] = result.pop("preservation")
    result["preservation"] = {"v1_raster_unchanged": True, "v1_manifest_byte_identical": True,
        "input_checksums_unchanged": True, "candidate_outputs_unchanged": True,
        "shared_manifest_advanced_from": 1, "shared_manifest_advanced_to": 2,
        "v1_manifest_archive": manifest["parent_stack"]["manifest_file"]}
    result["readiness"]["promotion"] = "approved by owner; published without retuning"
    result["promotion"] = "Performed after owner approval; exact COG copied, original v1 manifest archived, shared manifest advanced last"
    return result


def accepted_markdown(report: dict) -> bytes:
    # Reuse all scientific sections without changing the frozen builder or
    # retaining its now-inapplicable candidate-only/pending-review headings.
    original = diagnostics.markdown_report(report)
    body = original.split("## Approved configuration", 1)[1].rsplit("## Review and next step", 1)[0]
    publication = report["publication"]
    lines = ["# Accepted stack v2 scoring report", "", "Status: **accepted_owner_review** — published with all disclosed advisory limitations and review flags.", "",
        f"Raster: `{report['raster_file']}`", f"Raster SHA-256: `{report['raster_sha256']}`",
        f"Published manifest SHA-256: `{report['manifest_sha256']}`",
        f"Reviewed candidate manifest SHA-256: `{publication['reviewed_candidate']['manifest_sha256']}`", "",
        "## Owner decision", "", f"{report['owner_approval']['date']} — {report['owner_approval']['authority']}:", "",
        f"> {report['owner_approval']['statement']}", "",
        "The owner approved finalization after the diagnostic flags were disclosed. No thresholds were retuned. Owner acceptance is not independent field/agronomic validation, and the agent does not claim a visual inspection.", "",
        "## Post-copy technical acceptance", "",
        f"- {report['technical']['band_count']} bands; every full-resolution pixel and overview passed read-only validation of the copied COG.",
        f"- {report['technical']['scalar_checks']['compared_cells']} scalar reports match bitwise; the 38 v1 bands/records/tags are unchanged.",
        "- The published raster is byte-identical to the approved candidate. Candidate reproducibility evidence and all scientific findings remain unchanged.",
        f"- The original v1 manifest is archived byte-for-byte at `{report['preservation']['v1_manifest_archive']}`.", "",
        "## Approved configuration" + body,
        "## Phase 3 handoff", "",
        "Use the accepted shared manifest and common_stack_v2.tif. Resolve bands by name and restore point scoring with config.stack_v2.frozen_engine(manifest). Scores remain 0–1; missing JSON values are null; all values are nominal 1 km advisory estimates.", "",
        "Django/API/UI, PostGIS/PMTiles, roads/access infrastructure and optional external agronomic comparisons remain outside this Phase 2 finalization.", ""]
    return "\n".join(lines).encode("utf-8")


def finalize(args) -> dict:
    candidate_dir, work_dir = args.candidate_dir.resolve(), args.work_dir.absolute()
    processed = ROOT / "data/processed"
    if (".." in work_dir.parts or not work_dir.is_relative_to(ROOT / "data/interim") or work_dir == ROOT / "data/interim"
            or work_dir.is_relative_to(candidate_dir) or candidate_dir.is_relative_to(work_dir)):
        raise ValueError("publication must use a separate descriptive directory under data/interim")
    manifest, report, v1, approval = candidate_package(candidate_dir, args.approval)
    native_before = provenance.input_inventory(ROOT)
    scripts = publisher_files(ROOT)
    if not scripts or any(path.is_symlink() for path in scripts):
        raise ValueError("publisher source/snapshot files must be regular project files")
    promoted = published_manifest(manifest, approval, args.approval, candidate_dir, work_dir, processed, scripts)
    published_bytes = provenance.json_bytes(promoted)
    v1_bytes = (ROOT / manifest["parent_stack"]["manifest_file"]).read_bytes()
    shared = processed / "band_manifest.json"
    if shared.read_bytes() not in (v1_bytes, published_bytes):
        raise ValueError("shared manifest is neither the immutable v1 nor this approved v2; do not overwrite it")
    archive = processed / "band_manifest_v1.json"
    if archive.exists() and archive.read_bytes() != v1_bytes:
        raise ValueError("existing versioned v1 manifest differs; do not overwrite it")
    identity = {"candidate_manifest_sha256": provenance.sha256(candidate_dir / "band_manifest.json"),
                "candidate_raster_sha256": manifest["raster_sha256"], "approval_sha256": provenance.sha256(args.approval),
                "publisher_inputs": provenance.fingerprints(scripts, ROOT)}
    workspace = PublicationWorkspace(work_dir, processed, identity, scripts)
    protected = {path.resolve() for path in (args.approval, *(ROOT / name for name in manifest["build_inputs"]),
                 *(path for path in candidate_dir.rglob("*") if path.is_file()), *scripts)}
    all_targets = (*workspace.targets, workspace.state_path, *(path.with_name(path.name + ".tmp") for path in (*workspace.targets, workspace.state_path)))
    if any(path.resolve() in protected for path in all_targets):
        raise ValueError("publication output overlaps a protected candidate, approval or input")
    if workspace.state["status"] == "accepted":
        return {**validate_only(processed), "status": "already_finalized"}
    work_dir.mkdir(parents=True, exist_ok=True)
    workspace.checkpoint()
    log = []

    def note(message):
        print(message, flush=True)
        log.append(message)
        # The log is the sole intentionally overwritten, owned progress output.
        guard(work_dir / "finalize.log", ROOT)
        guard(work_dir / "finalize.log.tmp", ROOT)
        provenance.write_bytes(work_dir / "finalize.log", ("\n".join(log) + "\n").encode())
        workspace.checkpoint(work_dir / "finalize.log")

    for path in scripts:
        workspace.write(work_dir / "script_snapshot" / provenance.label(path, ROOT), path)
    note("Exact candidate, completed evidence, owner approval and all disclosed flags verified.")
    workspace.write(archive, v1_bytes)
    note("Original v1 manifest archived byte-for-byte before any shared-manifest advancement.")
    workspace.write(processed / "common_stack_v2.tif", ROOT / manifest["raster_file"])
    note("Approved COG copied without rewriting raster data; checking every copied pixel/overview and scalar report read-only.")
    acceptance = full_validation(promoted, v1, processed)
    accepted = accepted_report(report, promoted, approval, acceptance)
    payloads = {"band_manifest_v2.json": published_bytes, "stack_v2_validation.json": provenance.json_bytes(acceptance),
                "score_build_report.json": provenance.json_bytes(accepted), "score_build_report.md": accepted_markdown(accepted)}
    for name, payload in payloads.items():
        workspace.write(work_dir / name, payload)
    build_scores.verify_fingerprints(manifest["build_inputs"], ROOT, "publication input preservation")
    if provenance.input_inventory(ROOT) != native_before or fingerprint(args.approval)["sha256"] != identity["approval_sha256"]:
        raise ValueError("an input inventory entry or owner approval changed during publication")
    provenance.check_workspace(candidate_dir, ROOT, code_files=build_scores.code_files(ROOT))
    for path in scripts:
        if fingerprint(path) != identity["publisher_inputs"][provenance.label(path, ROOT)]:
            raise ValueError("publisher changed during finalization")
    for name in payloads:
        workspace.write(processed / name, work_dir / name)
    workspace.write(processed / "stack_v2_approval.json", args.approval)
    receipt_files = [processed / name for name in DELIVERABLES if name not in {"stack_v2_promotion.json", "band_manifest.json"}]
    receipt = {"status": "accepted", "stack_version": 2, "publisher": PUBLISHER, "owner_approval": approval,
        "reviewed_candidate": promoted["publication"]["reviewed_candidate"],
        "published_manifest_sha256": acceptance["manifest_sha256"], "raster_sha256": manifest["raster_sha256"],
        "published_files": provenance.fingerprints(receipt_files, ROOT),
        "input_checksums_unchanged": True, "native_input_inventory_unchanged": True, "candidate_outputs_unchanged": True,
        "v1_raster_unchanged": True, "v1_manifest_byte_identical": True,
        "publication_method": "Byte-exact COG copy, fresh full read-only acceptance, immutable v1 archive, owner approval; shared manifest advanced last",
        "scope": SCOPE}
    workspace.write(work_dir / "stack_v2_promotion.json", provenance.json_bytes(receipt))
    workspace.write(processed / "stack_v2_promotion.json", work_dir / "stack_v2_promotion.json")
    if shared.read_bytes() not in (v1_bytes, published_bytes) or archive.read_bytes() != v1_bytes:
        raise ValueError("shared/v1 archived manifest changed before activation; do not overwrite it")
    workspace.write(shared, published_bytes, replace_v1=True)
    stack_v2.assert_stack_matches_manifest(processed / "common_stack_v2.tif", promoted, v1_manifest=v1)
    note(f"Published {acceptance['band_count']} bands; {acceptance['scalar_checks']['compared_cells']} scalar checks matched. All scientific flags retained; no Phase 3 work performed.")
    workspace.state["status"] = "accepted"
    workspace.checkpoint()
    return {"status": "finalized", "band_count": acceptance["band_count"],
            "scalar_checks": acceptance["scalar_checks"]["compared_cells"], "raster_sha256": manifest["raster_sha256"],
            "manifest_sha256": acceptance["manifest_sha256"], "review_flags": approval["accepted_review_flags"]}


def validate_only(processed: Path) -> dict:
    manifest_path = processed / "band_manifest.json"
    manifest = provenance.load_json(manifest_path)
    stack_v2.validate_manifest(manifest)
    publication = manifest.get("publication", {})
    if manifest.get("status") != "accepted" or publication.get("approved") is not True:
        raise ValueError("read-only publication validation requires an owner-approved accepted v2")
    approval_path = ROOT / publication["approval_file"]
    if provenance.sha256(approval_path) != publication["approval_sha256"]:
        raise ValueError("published owner approval no longer matches its manifest")
    receipt = provenance.load_json(ROOT / publication["receipt_file"])
    if (receipt.get("status") != "accepted" or receipt.get("published_manifest_sha256") != provenance.sha256(manifest_path)
            or receipt.get("raster_sha256") != manifest["raster_sha256"]):
        raise ValueError("publication receipt does not identify this exact accepted package")
    build_scores.verify_fingerprints(receipt["published_files"], ROOT, "published deliverables")
    if manifest_path.read_bytes() != (processed / "band_manifest_v2.json").read_bytes():
        raise ValueError("shared and versioned v2 manifests differ")
    build_scores.verify_fingerprints(manifest["build_inputs"], ROOT, "published input provenance")
    build_scores.verify_fingerprints(publication["publisher_inputs"], ROOT, "publisher provenance")
    v1_path = ROOT / manifest["parent_stack"]["manifest_file"]
    if provenance.sha256(v1_path) != manifest["parent_stack"]["manifest_sha256"]:
        raise ValueError("versioned v1 archive no longer matches the immutable parent")
    v1, _ = build_scores.accepted_v1(ROOT / manifest["parent_stack"]["raster_file"], v1_path,
                                   ROOT / manifest["parent_acceptance"]["file"], ROOT)
    report = provenance.load_json(processed / "score_build_report.json")
    approval = provenance.load_json(approval_path)
    reviewed = publication["reviewed_candidate"]
    check_approval(approval, {"raster_sha256": reviewed["raster_sha256"]}, reviewed["manifest_sha256"], report)
    if publication["accepted_review_flags"] != approval["accepted_review_flags"] or reviewed["raster_sha256"] != manifest["raster_sha256"]:
        raise ValueError("published raster/review flags differ from owner approval")
    candidate_dir = (ROOT / reviewed["manifest_file"]).parent
    candidate, original, _, _ = candidate_package(candidate_dir, approval_path)
    expected = published_manifest(candidate, approval, approval_path, candidate_dir,
        ROOT / publication["workspace"], processed, publisher_files(ROOT))
    if provenance.json_bytes(expected) != manifest_path.read_bytes() or report["diagnostics"] != original["diagnostics"]:
        raise ValueError("publication changed approved scoring metadata or scientific evidence")
    if (report.get("status") != "accepted_owner_review" or report.get("manifest_sha256") != provenance.sha256(manifest_path)
            or report.get("owner_approval") != approval or receipt.get("owner_approval") != approval):
        raise ValueError("accepted report/receipt/approval are inconsistent")
    acceptance = full_validation(manifest, v1, processed)
    build_scores.verify_fingerprints(receipt["published_files"], ROOT, "read-only output preservation")
    return {"status": "passed_read_only", "band_count": acceptance["band_count"],
            "scalar_checks": acceptance["scalar_checks"]["compared_cells"], "raster_sha256": manifest["raster_sha256"],
            "manifest_sha256": provenance.sha256(manifest_path), "v1_preserved": True,
            "owner_review": "approved_with_disclosed_flags", "review_flags": publication["accepted_review_flags"],
            "scope": SCOPE}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-dir", type=Path, default=ROOT / "data/interim/score_build")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "data/interim/score_promote")
    parser.add_argument("--approval", type=Path, help="Explicit owner approval JSON naming both exact candidate hashes and all review flags")
    parser.add_argument("--validate-only", action="store_true", help="Revalidate the published package read-only; never rewrite files")
    args = parser.parse_args()
    if not args.validate_only and args.approval is None:
        parser.error("--approval is required for finalization; no implicit approval is inferred")
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetCacheMax(build_scores.GDAL_CACHE_BYTES)
    try:
        result = validate_only(ROOT / "data/processed") if args.validate_only else finalize(args)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Phase 2 finalization failed: {exc}\nNo automatic rollback or raster regeneration was performed; inspect publication_state.json before retrying.\n")
    print(provenance.json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
