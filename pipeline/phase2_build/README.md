# Offline Phase 2 scoring

Build a **candidate**, using the accepted v1 stack and owner-approved CSV/YAML
inputs. No Django, API, UI, PostGIS, PMTiles, native-raster regeneration or
publication is performed.

From the repository root, with the approved interpreter:

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python pipeline/phase2_build/build_scores.py \
  --crop-requirements data/crops/crop_requirments.csv
```

Paths can be overridden with the command's `--help` options. The misspelled crop
filename is the explicitly approved source, not a silently renamed copy.
Candidate output directories must be descriptive children of `data/interim/`.
The builder reruns its own outputs but refuses unfamiliar or externally edited
deliverables. It never edits the original CSVs, grid, stack v1 or shared manifest.

## Outputs

Under `data/interim/score_build/`:

- `common_stack_v2.tif`: candidate COG; 60 bands for all six active crops.
- `common_stack_v2_source.tif`: intermediate assembly with verified overviews.
- `band_manifest.json`: candidate v2 metadata, frozen rules/units/calibration,
  factor lookups, original crop IDs and input/code SHA-256 fingerprints.
- `band_manifest_v1.json`: original v1 manifest, archived byte-for-byte.
- `technical_validation.json`: every-pixel, every-overview and scalar evidence.
  It is checkpointed before diagnostics with an explicit incomplete status, then
  finalized only after reproducibility and input-preservation checks pass.
- `score_build_report.md` / `.json`: distributions, missing inputs, gates,
  geographic checks, crop-pair Spearman correlations, limiting shares,
  confidence-edge effects and approved class-change sensitivities.
- `reproducibility/`: independently rebuilt TIFFs; final COG and regenerated
  manifest must be byte-identical in the recorded runtime.
- `script_snapshot/`, `build.log`, `build_state.json`: exact code/grid snapshots,
  step log and narrow output ownership. Snapshots are evidence; run the command
  from `pipeline/phase2_build/`, not from the snapshot directory.

Recheck an existing candidate without rewriting files:

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python pipeline/phase2_build/build_scores.py \
  --validate-only
```

## Acceptance and interpretation

The first 38 v1 band records, full-resolution float32 bits and per-band provenance
tags remain unchanged. All appended scores are 0–1, with physical `-9999` nodata.
Categorical IDs use mode overviews; continuous scores use average overviews.
Every overview level is regenerated directly from full resolution.

The command budgets 1 GiB for GDAL's decoded-block cache. After every assembly
pixel is bit-verified against the COG, statistics/diagnostics read its faster
band-interleaved intermediate. This does not alter the candidate or weaken COG,
overview or scalar checks. Interrupted steps can be rerun with the same command;
only previously owned outputs are overwritten, and incomplete evidence is not
treated as a finished or approved candidate.

Single-cell reports use `config.stack_v2.frozen_engine(manifest)`, so subsequent
CSV edits cannot silently change the explanation of an old stack. Numerical
acceptance recomputes every scored pixel in differently sized windows and checks
32 seeded random inside-country cells plus available edge/exclusion cases.

Percent rainfall reliability is the approved Normal approximation for rainfall
reaching a configured fraction of seasonal **mean**, not a probability of meeting
crop water demand. Phosphorus remains mg/100kg; EthioSIS nutrients remain mg/kg.
Chickpea uses its approved Kiremt recharge/residual-moisture proxy. Unsourced
coarse fragments and coverage-radius sensitivity are disclosed and skipped.

Review flags are reported, **never fixed by retuning thresholds automatically**.
Coverage-edge diagnostics separate ordinary input variation from the bounded
confidence-blend effect; they do not claim visual or agronomic validation.
Original excluded/gapped capital cells remain visible beside disclosed nearby
valid diagnostic cells selected by grid distance, not by crop ranking.

Promotion requires owner review of the exact candidate raster/manifest hashes
and the observed diagnostics. Before any future shared-manifest advancement,
preserve a byte-identical versioned v1 manifest in `data/processed/`. This command
intentionally has no promotion switch.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m unittest discover \
  -s pipeline/phase2_build/tests -t .
```

Tests use invented data under `/tmp/opencode` and never regenerate real inputs.
