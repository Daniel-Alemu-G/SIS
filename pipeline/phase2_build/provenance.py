"""Checksums and narrow output ownership; never overwrite unfamiliar work."""

from __future__ import annotations

import hashlib
import json
import math
import shutil
from pathlib import Path

BUILD_SCRIPT = "pipeline/phase2_build/build_scores.py"
STATE_FILE = "build_state.json"
OUTPUT_FILES = (
    "common_stack_v2_source.tif", "common_stack_v2.tif", "common_stack_v2.tmp.tif",
    "band_manifest_v1.json", "band_manifest.json", "technical_validation.json",
    "score_build_report.json", "score_build_report.md", "build.log",
    "reproducibility/common_stack_v2_source.tif", "reproducibility/common_stack_v2.tif",
    "reproducibility/common_stack_v2.tmp.tif",
)


def sha256(path: Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def label(path: Path, root: Path) -> str:
    path, root = Path(path).resolve(), Path(root).resolve()
    return path.relative_to(root).as_posix() if path.is_relative_to(root) else str(path)


def fingerprints(paths, root: Path) -> dict:
    return {label(path, root): {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted({Path(path).resolve() for path in paths})}


def input_inventory(root: Path) -> dict:
    """Unused native inputs are inventoried, not read/reprocessed for scoring."""
    return {label(path, root): {"size_bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for path in sorted((root / "data/inputs").rglob("*")) if path.is_file()}


def json_bytes(value) -> bytes:
    # No timestamps, hostnames or file mtimes in reproducible artifacts.
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def load_json(path: Path) -> dict:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"{path}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    def constant(value):
        raise ValueError(f"{path}: nonfinite JSON value {value!r}")

    def number(value):
        result = float(value)
        if not math.isfinite(result):
            constant(value)
        return result

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_constant=constant, parse_float=number)


def write_bytes(path: Path, data: bytes) -> None:
    temporary = path.with_name(path.name + ".tmp")
    if path.is_symlink() or temporary.is_symlink():
        raise ValueError(f"refusing a symlink output: {path}")
    temporary.write_bytes(data)
    temporary.replace(path)


def write_json(path: Path, value) -> None:
    write_bytes(path, json_bytes(value))


def archive_manifest(source: Path, destination: Path) -> None:
    """Archive exact bytes, not reserialized JSON; an archive is immutable."""
    if source.resolve() == destination.resolve() or destination.is_symlink():
        raise ValueError("v1 manifest cannot be archived in place or through a symlink")
    original = source.read_bytes()
    if destination.exists():
        if not destination.is_file() or destination.read_bytes() != original:
            raise ValueError("existing v1 manifest archive differs; do not overwrite it")
        return
    write_bytes(destination, original)
    if destination.read_bytes() != original:
        raise ValueError("v1 manifest archive is not byte-identical")


def workspace_targets(managed) -> tuple[str, ...]:
    """Every owned write target, including atomic-write scratch files."""
    managed = tuple(managed)
    return (*managed, STATE_FILE, *(name + ".tmp" for name in (*managed, STATE_FILE)))


def check_workspace(work_dir: Path, root: Path, *, code_files=()) -> tuple[Path, tuple[str, ...]]:
    """Validate every target before claiming a fresh or previously owned step."""
    work_dir, root = work_dir.resolve(), root.resolve()
    interim = root / "data/interim"
    if not work_dir.is_relative_to(interim) or work_dir == interim:
        raise ValueError("candidate outputs must use a descriptive directory under data/interim")
    code_files = tuple(Path(path).resolve() for path in code_files)
    if any(not path.is_relative_to(root) for path in code_files):
        raise ValueError("script snapshots must come from files within the project")
    snapshots = tuple("script_snapshot/" + label(path, root) for path in code_files)
    if any(Path(name).is_absolute() or ".." in Path(name).parts for name in snapshots):
        raise ValueError("script snapshots must come from files within the project")
    managed = (*OUTPUT_FILES, *snapshots)
    # Include atomic-write scratch names in the safety checks, but do not claim
    # arbitrary files or clean unrelated contents of this existing directory.
    checked = workspace_targets(managed)
    for name in checked:
        path = work_dir / name
        parents = path.parents[:len(path.relative_to(work_dir).parts) - 1]
        if path.is_symlink() or any(parent.is_symlink() for parent in parents) or not path.resolve().is_relative_to(work_dir):
            raise ValueError(f"refusing a symlink/escaped output: {path}")
        if path.exists() and not path.is_file():
            raise ValueError(f"output path is not a file: {path}")
    state_path = work_dir / STATE_FILE
    if not state_path.exists():
        if any((work_dir / name).exists() for name in checked if name != STATE_FILE):
            raise ValueError("existing candidate outputs have no known owner; do not overwrite them")
    else:
        state = load_json(state_path)
        if not isinstance(state, dict) or state.get("build_script") != BUILD_SCRIPT or type(state.get("stack_version")) is not int or state["stack_version"] != 2 or state.get("candidate_only") is not True:
            raise ValueError("existing workspace belongs to an unfamiliar builder")
        previous = state.get("managed_files")
        if not isinstance(previous, list) or not set(OUTPUT_FILES) <= set(previous) or any(name not in managed for name in previous):
            raise ValueError("workspace output ownership differs; inspect it before rerunning")
        saved = state.get("completed_outputs")
        if not isinstance(saved, dict) or any(name not in previous for name in saved):
            raise ValueError("workspace completed-output fingerprints are invalid")
        for name, expected in saved.items():
            path = work_dir / name
            if not path.is_file() or {"size_bytes": path.stat().st_size, "sha256": sha256(path)} != expected:
                raise ValueError(f"previous candidate output changed outside this step: {path}")
        for name in set(managed) - set(previous):
            if (work_dir / name).exists():
                raise ValueError(f"newly claimed output already exists: {work_dir / name}")
    return work_dir, managed


def begin_workspace(work_dir: Path, managed: tuple[str, ...]) -> None:
    work_dir.mkdir(parents=True, exist_ok=True)
    write_json(work_dir / STATE_FILE, {"build_script": BUILD_SCRIPT, "stack_version": 2, "candidate_only": True,
                                     "managed_files": list(managed), "completed_outputs": {}, "status": "building"})


def finish_workspace(work_dir: Path, managed: tuple[str, ...]) -> None:
    completed = {name: {"size_bytes": (work_dir / name).stat().st_size, "sha256": sha256(work_dir / name)}
                 for name in managed if (work_dir / name).is_file()}
    write_json(work_dir / STATE_FILE, {"build_script": BUILD_SCRIPT, "stack_version": 2, "candidate_only": True,
                                     "managed_files": list(managed), "completed_outputs": completed,
                                     "status": "awaiting_owner_review"})


def snapshot_code(code_files, work_dir: Path, root: Path) -> None:
    for path in code_files:
        destination = work_dir / "script_snapshot" / label(path, root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        if sha256(path) != sha256(destination):
            raise ValueError(f"script snapshot checksum differs: {destination}")
