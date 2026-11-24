# Phase 1 progress

## Done

- `P1-01` — inventoried all inputs and recorded approved policies plus unresolved phosphorus provenance.
- `P1-03` — Python 3.12 `.venv` passes the Phase 1 import smoke test.
- `P1-04` — canonical-grid validator passes 6 unit tests.
- `P1-05` — synthetic float32/`-9999` COG preserves average and mode overviews.
- `P1-06` — mask (1,338,026 cells), 1,148-row woreda lookup, and stable ID raster pass rerun checks; two tiny woredas contain no cell centre.
- `P1-07` — all 12 EthioSIS bands window-copy to 1812 columns with no interpolation; per-nutrient null audits pass.
- `P1-08` — converted 269,774,816 valid 90 m percentage pixels to degrees before averaging; all 1,338,026 masked cells are valid and span 0–49.62°.

## In progress

- `P1-09` — stage exact-grid SoilGrids bands except gated phosphorus, plus `elevation_m`.

## Decisions

- Grid: use only `pipeline/grid.yaml` via `config.grid`—EPSG:4326, pixel-is-area, 1/120°, 1812×1380.
- Stack: append-only `common_stack_v{N}.tif`; access bands by lowercase snake_case name through `band_manifest.json`, never by index.
- Storage: one float32 COG with `-9999` physical nodata; manifest retains logical types. Continuous overviews use average; categorical/ID/mask use mode.
- Masking: all bands are nodata outside `ethiopia_mask`; the mask itself remains meaningful 0/1. `woreda_id` uses 0 for none inside Ethiopia.
- Units: soil texture %, CEC cmol(c)/kg, SOC g/kg, `sg_n` mg/kg, CFVO vol %, EthioSIS nutrients mg/kg by source header, elevation m, slope degrees, rainfall mm/CV fraction, temperature °C, distance km. Phosphorus remains gated on provenance/unit/scaling.
- Nitrogen: `sg_n` is the sole application/scoring input and is never rescaled again; do not create EthioSIS N, P, or OC bands.
- EthioSIS: window-copy the first 1812 columns without interpolation; preserve nutrient-specific nulls and use confidence only for sample proximity.
- Slope: convert each valid 90 m percentage pixel with `degrees(arctan(percent/100))` before averaging; keep both 90 m terrain rasters separate for exact clicks.
- Climate: use kiremt Jun–Sep, belg Feb–May, bega Oct–Jan; reuse the nine-band seasonal rainfall product rather than monthlies.
- Stores: per-cell numbers go in the COG; named woredas, roads, capitals, and samples go in PostGIS SRID 4326 with spatial indexes.
- Boundaries/access: use approved v04 polygons; show regions/woredas only. Core v1 has 38 bands; all three access-distance bands are deferred.
- Roads/centres: use the approved motor-road allowlists, preserve source names/IDs, and label nearest-capital context as derived. Capitals are incomplete centre coverage, not towns or markets.
- Runtime/delivery: use a Python 3.12 project venv. Named-object loaders and PMTiles wait for a test database and compatible tooling; no untestable assets.

## Blockers

- `P1-02`: owner must provide the exact phosphorus product, unit, scaling, and back-transformation before `sg_p` or stack v1 can be finalized.
- `P1-14`: PostGIS and PMTiles work remains deferred until a test database and compatible tooling are supplied.
