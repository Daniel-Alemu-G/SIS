# Finalize the reviewed Phase 2 candidate

This command **copies** an already built COG. It does not regenerate any raster,
retune a threshold, change the scoring implementation or start Phase 3.

An explicit owner approval JSON must name the exact candidate raster and manifest
SHA-256 values and every disclosed diagnostic flag. Approval is never a command
default. Keep the decision under `data/interim/score_promote/owner_approval.json`:

```json
{
  "approved": true,
  "authority": "Project owner",
  "date": "YYYY-MM-DD",
  "statement": "The owner's actual approval statement",
  "scope": "Phase 2 scoring stack v2 finalization only",
  "candidate_raster_sha256": "<exact reviewed raster SHA-256>",
  "candidate_manifest_sha256": "<exact reviewed manifest SHA-256>",
  "accepted_review_flags": ["<every report warning code, in report order>"]
}
```

With the approved project interpreter, from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python pipeline/phase2_promote/finalize_scores.py \
  --approval data/interim/score_promote/owner_approval.json

PYTHONDONTWRITEBYTECODE=1 .venv/bin/python pipeline/phase2_promote/finalize_scores.py \
  --validate-only
```

## Publication contract

- Verify the completed candidate, review report, approval and consumed-input/code
  fingerprints before claiming outputs. Retain all six-crop scientific findings.
- Archive the original v1 manifest **byte-for-byte** at
  `data/processed/band_manifest_v1.json` before the shared manifest advances.
- Copy `common_stack_v2.tif` without changing its approved bytes. Recheck every
  full-resolution pixel, overview, v1 record/tag and scalar report read-only.
- Stage the accepted metadata/reports and exact publisher snapshots under
  `data/interim/score_promote/`, then copy the approved package to `data/processed/`.
- Publish `band_manifest_v2.json`, `stack_v2_approval.json`,
  `stack_v2_validation.json`, `stack_v2_promotion.json`, and
  `score_build_report.md`/`.json`. Advance `band_manifest.json` **last**; it is
  byte-identical to the versioned v2 manifest.
- The published manifest has a new metadata hash because status, publication
  paths, the stable v1 archive reference and owner approval are recorded. Its
  approval record still identifies the **original reviewed candidate hash**.
  Band records, frozen scoring/calibration, input provenance and COG bytes are
  unchanged. The original candidate and reports remain untouched.
- Refuse unfamiliar, externally modified, aliased or conflicting outputs. Narrow
  checkpoints allow safe interrupted runs; a completed rerun only revalidates.
- `--validate-only` changes no files and verifies the package, original candidate,
  owner decision, retained diagnostics and full scoring/storage contract.

Owner acceptance is advisory acceptance, not proof of external field/agronomic
validation. Do not hide flags, assert a visual inspection that did not occur, or
silently replace owner thresholds.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m unittest discover \
  -s pipeline/phase2_promote/tests -t .
```

Tests use synthetic files under `/tmp/opencode`, never real project outputs.
