# Phase 1 input inventory report

**Audit date:** 2026-11-23

**Status:** complete read-only inventory; raster-foundation decisions approved; AfricanSoilGrids phosphorus approved as-is in mg/100kg; raster scientific gates closed

**Phosphorus follow-ups:** 2026-11-24 (initial unit confirmation in section 13; final owner approval in section 14)

**Canonical grid source:** `pipeline/grid.yaml`

## 1. Scope and method

- Audited all 219 files under `data/inputs/` (2,442,948,648 bytes): 198 GeoTIFFs, 7 GeoJSONs, 2 CSVs, 1 GeoPackage, 7 PAM XML sidecars, 2 JSON files, 1 YAML file and 1 text file.
- All inputs were opened read-only with GDAL/OGR or SQLite read-only mode. Size and nanosecond modification-time snapshots were identical before and after both completed audits.
- Raster min/max/counts are full scans. Medians are exact for common-grid rasters and the integer DEM; the two large float32 90 m medians use a systematic 1/256 sample (every 16th row and column), while their min/max/counts are full scans.
- Common-grid inside-country statistics use an in-memory cell-centre rasterization of `eth_admin0.geojson`; the mask contains 1,338,026 cells. This audit created no project data derivative.
- The canonical YAML was read directly for this audit because the current system interpreter lacks `affine`, so `config/grid.py` cannot import. This is an environment finding, not a second grid definition.
- Input GeoTIFF band numbers below are inventory-only. Production access remains by band description/name.

## 2. Executive findings

1. All intended 1 km inputs match the canonical CRS, affine transform, pixel-is-area convention and 1812 × 1380 shape **except** `ethiosis_IDW_estimates.tif`, which is 1813 × 1380. Its extra eastern column has NaN in all 11 nutrient bands and zero confidence.
2. Both slope rasters contain **percentage slope**, confirmed by source code (`slopeFormat="percent"`) and observed ranges. The approved output remains mean degrees: convert valid 90 m pixels with `degrees(arctan(percent / 100))` before common-grid averaging.
3. `sg_n` is already in mg/kg and is the sole application/scoring nitrogen source. Its inside-country median is 1,247.783 mg/kg. No additional scaling or EthioSIS nitrogen input is permitted.
4. The IDW raster contains `ca, cu, fe, mg, mn, mo, k, si, na, s, zn, confidence`; it has no OC, N or P band. Positive confidence does not guarantee that every individual nutrient is present (84,204 positive-confidence cells lack molybdenum).
5. The seasonal rainfall COG and sidecar agree, all CV bands equal `std / mean`, and all 181 monthly files form a gap-free 2011-09 through 2026-09 sequence. Reuse the seasonal product; do not rebuild it.
6. Temperature values near -9929 occur only outside the canonical Ethiopia mask. All three bands have zero values outside the -20 to 70 °C sanity range inside Ethiopia, so the inventory gives no reason to regenerate temperature.
7. The roads GeoPackage has one mixed-geometry `roads` layer (450,047 records: 447,034 lines, 2,617 points, 396 polygons), no places layer, no `fclass` or `ref`, and typed string IDs. Only 10,290 records have any preserved name field populated.
8. `eth_admincapitals.geojson` has 857 valid points and all names populated. `loc_pcode` is empty throughout; 15 level/P-code keys repeat and the national Addis Ababa record has no level-matched P-code. A persisted source-aware record ID is required.
9. The woreda GeoJSON has 1,148 valid polygons/multipolygons, no market fields, and no populated centre coordinates. The separate administrative gazetteer has 1,148 populated coordinate pairs plus two extra P-codes; neither file establishes market locations.
10. On 2026-11-23 the owner approved the storage/nodata, EthioSIS alignment/missingness, boundary, road-class, access-band, runtime and deferral decisions in section 11. The final v1 stack remains blocked only on the requested phosphorus provenance details; approved preparatory transformations and tests may proceed.
11. The initial 2026-11-24 phosphorus unit confirmation established **mg/100kg** and a zero-mismatch interim masked copy. At that point product/transform approval remained open; the subsequent owner confirmation in finding 12 supersedes that initial gate.
12. The owner subsequently identified **AfricanSoilGrids**, explicitly distinguished it from global SoilGrids, and approved the supplied phosphorus values **as-is, with no scaling or change**. `sg_p` retains mg/100kg, source gaps and all valid values; no rescaling or back-transformation is applied. This closes `P1-02` and enables final v1 assembly.

## 3. Canonical grid and environment

| Item | Observed value |
|---|---|
| CRS | EPSG:4326 |
| Pixel convention | pixel-is-area (`AREA_OR_POINT=Area` on all common-grid rasters) |
| Shape | 1812 columns × 1380 rows |
| Affine (GDAL order) | `(32.9, 1/120, 0, 14.9, 0, -1/120)` |
| Extent | west 32.9, east 48.0, south 3.4, north 14.9 |
| YAML nodata | `-9999` |
| System interpreter | `/usr/bin/python3` 3.12.3 |
| Available Python modules | GDAL/`osgeo`, numpy, pyproj, PyYAML, scipy |
| Missing required/proposed modules | affine, rasterio, geopandas, shapely, Django, psycopg |
| Available CLIs | `gdalinfo`, `ogrinfo`, `ogr2ogr`, Node 24.18.0, npm |
| Missing CLIs | `psql`, `pg_isready`, `tippecanoe`, `pmtiles` |
| Scratch space observed | about 61.9 GB free under `/tmp/opencode` |

Do not install packages into system Python. Confirm an isolated interpreter and dependency manifest first. No database target was available to audit.

## 4. Raster geometry and storage inventory

| Dataset | Shape (rows × cols) | CRS/resolution | Grid | Physical storage | Nodata | Bands/overviews |
|---|---:|---|---|---|---|---|
| `data/inputs/Elevation/ethiopia_Copernicus-DEM_1KM.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 1812×1 | -9999.0 | 1; overview levels [0] |
| `data/inputs/Elevation/ethiopia_copernicus_dem_90m.tif` | 14,400 × 18,000 | EPSG:4326, 3 arc-sec (~90 m) | separate by design | Float32, DEFLATE, block 2048×2048 | None | 1; overview levels [0] |
| `data/inputs/EthioSIS/ethiosis_IDW_estimates.tif` | 1,380 × 1,813 | EPSG:4326, 1/120° with extra east column | shape mismatch only | Float32, DEFLATE, block 256×256 | NaN | 12; overview levels [0] |
| `data/inputs/Landcover/WorldCover_Ethiopia_1km.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Byte, DEFLATE, block 1812×4 | 0.0 | 1; overview levels [0] |
| `data/inputs/Rainfall/seasonal/rainfall_seasonal_1km.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 512×512, COG declared | -9999.0 | 9; overview levels [2] |
| `data/inputs/Slope/slope_1km.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/Slope/slope_albers_90m.tif` | 14,905 × 18,527 | custom Albers, 90 m | separate by design | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/cec_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/cfvo_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/clay_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/n_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/p_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/phh2o_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/sand_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/silt_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/SoilGrids/soc_0-30cm_mean_1km_4326.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, DEFLATE, block 256×256 | -9999.0 | 1; overview levels [0] |
| `data/inputs/temprature/tmean_seasonal.tif` | 1,380 × 1,812 | EPSG:4326, 1/120° | exact | Float32, LZW, block 256×256 | -9999.0 | 3; overview levels [0] |

Only the seasonal rainfall file declares COG layout and has internal overviews (906 × 690 and 453 × 345). Other inputs may be tiled/striped but are not accepted as serving COGs without validation/conversion. Input unit tags, scales and offsets are empty on every inspected raster; units therefore depend on documented provenance and CSV headers, not embedded metadata.

### 4.1 Common-grid values inside Ethiopia

| File/band | Unit basis | Min | Median | Max | Missing mask cells | Finding |
|---|---|---:|---:|---:|---:|---|
| `cec_0-30cm_mean_1km_4326.tif` / `cec_0-30cm_mean_1km_4326` | cmol(c)/kg | 10.901 | 28.707 | 84.253 | 6,213 | ok |
| `cfvo_0-30cm_mean_1km_4326.tif` / `cfvo_0-30cm_mean_1km_4326` | vol % | 1.802 | 16.96 | 56.901 | 6,213 | ok |
| `clay_0-30cm_mean_1km_4326.tif` / `clay_0-30cm_mean_1km_4326` | % | 14.585 | 34.965 | 69.441 | 6,213 | ok |
| `n_0-30cm_mean_1km_4326.tif` / `n_0-30cm_mean_1km_4326` | mg/kg | 446.026 | 1,247.783 | 8,191.428 | 5,416 | confirmed mg/kg; values exceed rough 500–5000 range |
| `p_0-30cm_mean_1km_4326.tif` / unnamed singleton | mg/100kg (owner, 2026-11-24) | 174.626 | 1,596.827 | 9,091.994 | 5,416 | AfricanSoilGrids; owner approved values as-is; 2,200 cells above provisional range, not clipped |
| `phh2o_0-30cm_mean_1km_4326.tif` / `phh2o_0-30cm_mean_1km_4326` | pH | 4.931 | 7.121 | 8.86 | 6,213 | ok |
| `sand_0-30cm_mean_1km_4326.tif` / `sand_0-30cm_mean_1km_4326` | % | 3.36 | 33.886 | 69.007 | 6,213 | ok |
| `silt_0-30cm_mean_1km_4326.tif` / `silt_0-30cm_mean_1km_4326` | % | 10.22 | 30.861 | 65.287 | 6,213 | ok |
| `soc_0-30cm_mean_1km_4326.tif` / `soc_0-30cm_mean_1km_4326` | g/kg | 4.773 | 18.484 | 91.783 | 6,962 | ok |
| `ethiosis_IDW_estimates.tif` / `ca` | mg/kg (CSV/header basis) | 31.941 | 4,194.953 | 25,614.318 | 375,018 | 107 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `cu` | mg/kg (CSV/header basis) | 0.0025 | 2.67 | 12.907 | 374,958 | 52 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `fe` | mg/kg (CSV/header basis) | 0.0498279 | 94.536 | 556.164 | 374,977 | 66 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `mg` | mg/kg (CSV/header basis) | 5.885 | 574.24 | 3,478.871 | 374,958 | 47 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `mn` | mg/kg (CSV/header basis) | 1.513 | 129.138 | 540.683 | 375,064 | 153 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `mo` | mg/kg (CSV/header basis) | 0.000157 | 0.491302 | 939.1 | 453,326 | 84,204 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `k` | mg/kg (CSV/header basis) | 9.237 | 337.178 | 2,365.372 | 375,017 | 106 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `si` | mg/kg (CSV/header basis) | 0.1 | 476.993 | 1,235.832 | 375,799 | 888 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `na` | mg/kg (CSV/header basis) | 0.511855 | 37.032 | 5,620.801 | 375,051 | 140 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `s` | mg/kg (CSV/header basis) | 0.301474 | 11.879 | 7,259.35 | 375,106 | 195 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `zn` | mg/kg (CSV/header basis) | 0.003067 | 1.339 | 38.572 | 375,117 | 206 missing with positive confidence |
| `ethiosis_IDW_estimates.tif` / `confidence` | fraction | 0 | 0.646738 | 1 | 0 | ok |
| `ethiopia_Copernicus-DEM_1KM.tif` / `ethiopia_Copernicus-DEM_1KM` | m | -125 | 1,146.579 | 4,390.596 | 5 | ok |
| `slope_1km.tif` / `slope_1km` | % slope (native input) | 0 | 5.827 | 138.584 | 0 | percent, not degrees; do not use as slope_mean_deg |
| `WorldCover_Ethiopia_1km.tif` / `Map` | class code | 10 | 20 | 90 | 4,501 | legal observed codes 10,20,30,40,50,60,80,90 |
| `rainfall_seasonal_1km.tif` / `kiremt_mean_mm` | mm | 4.838 | 413.33 | 2,303.243 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `kiremt_std_mm` | mm | 1.178 | 89.978 | 298.116 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `kiremt_cv` | fraction | 0.0515514 | 0.228974 | 1.71 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `belg_mean_mm` | mm | 48.916 | 251.534 | 1,089.736 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `belg_std_mm` | mm | 28.84 | 91.964 | 251.547 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `belg_cv` | fraction | 0.0786987 | 0.362283 | 0.767042 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `bega_mean_mm` | mm | 17.521 | 143.636 | 786.183 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `bega_std_mm` | mm | 9.623 | 72.878 | 228.332 | 0 | ok |
| `rainfall_seasonal_1km.tif` / `bega_cv` | fraction | 0.118402 | 0.503594 | 1.128 | 0 | ok |
| `tmean_seasonal.tif` / `tmean_kiremt` | °C | 9.485 | 24.996 | 39.502 | 0 | 0 out-of-range inside; 21,078 outside-country outliers |
| `tmean_seasonal.tif` / `tmean_belg` | °C | 10.098 | 26.073 | 36.469 | 0 | 0 out-of-range inside; 21,078 outside-country outliers |
| `tmean_seasonal.tif` / `tmean_bega` | °C | 9.705 | 24.405 | 35.095 | 0 | 0 out-of-range inside; 21,078 outside-country outliers |

Land-cover counts over the raster valid footprint: `10`=184,331, `20`=512,514, `30`=272,682, `40`=252,634, `50`=1,980, `60`=97,871, `80`=9,813, `90`=9,151. Codes 70, 95 and 100 are absent; 0 is source nodata. There are 4,501 land-cover nodata cells inside the administrative mask.

### 4.2 90 m terrain rasters

| File | CRS/grid | Nodata | Valid count | Min | Median | Max | Notes |
|---|---|---|---:|---:|---:|---:|---|
| `data/inputs/Elevation/ethiopia_copernicus_dem_90m.tif` | EPSG:4326, 3 arc-sec | None | 259,200,000 | -154.818 | 548.886 | 4,511.379 | No declared nodata; 85,406,238 zeros exist and must not automatically be reclassified. |
| `data/inputs/Slope/slope_albers_90m.tif` | custom Albers, 90 m | -9999.0 | 269,774,816 | 0 | 0.846405 | 1,084.304 | Native percent slope; 6,370,119 pixels are -9999 nodata. Approximate median; max 1084.304%. |

### 4.3 EthioSIS alignment and coverage

- Shape is 1813 × 1380; CRS and all affine components match the canonical grid. The extra eastern column has 0 finite nutrient values; all 1,380 confidence values are zero.
- Confidence is present for every cell and ranges 0–1. Nutrient bands use NaN nodata. No nutrient is valid where confidence equals zero.
- Missing nutrient cells with positive confidence: `ca` 107, `cu` 52, `fe` 66, `mg` 47, `mn` 153, `mo` 84,204, `k` 106, `si` 888, `na` 140, `s` 195, `zn` 206.
- This disproves the current strict rule that positive confidence implies every nutrient exists. It does not justify imputation.

### 4.4 Rainfall source series

- 181 monthly CHIRPS files, one per month from September 2011 through September 2026, with no duplicate or missing month.
- Every monthly file is one float32 band, EPSG:4326, pixel-is-area, canonical 1812 × 1380 grid, -9999 nodata, DEFLATE/BAND storage, and no band description. Full-scan per-file statistics were inspected; across the series the global valid minimum is 0 mm and maximum is 1,026.897 mm.
- The existing 9-band seasonal COG has exact canonical geometry. Its companion JSON exactly matches all band descriptions. Kiremt and belg use 2012–2026 (15 years); bega uses 2011–2025 (15 years).
- Full-cell checks found zero `cv != std / mean` mismatches for all 2,323,440 valid cells in each season.
- Disposition: monthly files are upstream provenance only for Phase 1. Integrate `rainfall_seasonal_1km.tif`; do not reaggregate the monthlies.

## 5. Vector inventory

| File/layer | Count | Geometry | CRS | Validity/identity findings | Intended disposition |
|---|---:|---|---|---|---|
| `data/inputs/Admin boundaries/eth_admin0.geojson` / `eth_admin0` | 1 | POLYGON 1 | EPSG:4326 | 0 missing/empty; 0 invalid polygons | country mask and provenance |
| `data/inputs/Admin boundaries/eth_admin1_regions.geojson` / `eth_admin1` | 15 | POLYGON 13, MULTIPOLYGON 2 | EPSG:4326 | 0 missing/empty; 0 invalid polygons | region lookup/display source |
| `data/inputs/Admin boundaries/eth_admin2_zone.geojson` / `eth_admin2` | 107 | POLYGON 101, MULTIPOLYGON 6 | EPSG:4326 | 0 missing/empty; 0 invalid polygons | parent metadata only; no production zone display geometry |
| `data/inputs/Admin boundaries/eth_admin3_woreda.geojson` / `eth_admin3` | 1,148 | POLYGON 1,125, MULTIPOLYGON 23 | EPSG:4326 | 0 missing/empty; 0 invalid polygons; 1,148 unique adm3_pcode values; centre fields empty | woreda raster IDs and full-resolution lookup if owner confirms provenance/resolution |
| `data/inputs/Admin boundaries/eth_admincapitals.geojson` / `eth_admincapitals` | 857 | POINT 857 | EPSG:4326 | 0 missing/empty; 0 invalid polygons; loc_pcode empty; 15 repeated level/P-code keys; 1 missing matched P-code | approved centre lookup and derived unnamed-road context |
| `data/inputs/Admin boundaries/eth_adminlines.geojson` / `eth_adminlines` | 4,290 | LINESTRING 4,289, MULTILINESTRING 1 | EPSG:4326 | 0 missing/empty; 0 invalid polygons | ancillary; inventory only |
| `data/inputs/Admin boundaries/eth_adminpoints.geojson` / `eth_adminpoints` | 1,271 | POINT 1,271 | EPSG:4326 | 0 missing/empty; 0 invalid polygons | ancillary; inventory only |
| `data/inputs/Road data/eth_roads_osm_gpkg/roads.gpkg` / `roads` | 450,047 | LINESTRING 447,034, POINT 2,617, POLYGON 396 | EPSG:4326 | 0 missing/empty; IDs unique; 447,034 line candidates; no places/ref/fclass | normalize approved lines to PostGIS/PMTiles after class sign-off |

### 5.1 Administrative schemas and provenance

- All administrative geometries are EPSG:4326, version `v04`, `valid_on=2025-01-01`; inspected polygons are OGC-valid. The files contain no licence/source declaration sufficient to confirm publication rights or that lookup polygons are unsimplified/full-resolution.
- Region: 15 records (`adm1_name`, `adm1_pcode`, country parents and metadata). Zone: 107. Woreda: 1,148 with unique `adm3_pcode` values (maximum length 8), parent zone/region names and P-codes, and polygon area metadata.
- `eth_admin3_gzt.csv` has 1,150 unique P-code rows; 1,148 have coordinate pairs matching polygon P-codes, while `ET042199` and `ET050905` exist only in the gazetteer. These are administrative gazetteer points, not verified markets.
- `eth_admincapitals.geojson` levels: national 1, regional 12, zonal 84, woreda 760. All 857 names/geometries exist; `loc_pcode` and all country fields are empty. Distinct records must not be merged only because a level/P-code repeats.
- `eth_adminpoints.geojson` (1,271 points) and `eth_adminlines.geojson` (4,290 lines) are unrequested ancillary products. Do not silently substitute them for the approved capitals or polygon sources.

### 5.2 Roads

- Source metadata: OpenStreetMap/Geofabrik Ethiopia snapshot 2026-10-12, exported 2026-10-13, ODbL-1.0, EPSG:4326.
- Fields: `id,name,name_en,name_am,highway,surface,smoothness,width,lanes,oneway,bridge,layer,source,adm0_pcode,adm0_name,adm1_pcode,adm1_name,adm2_pcode,adm2_name,adm3_pcode,adm3_name,adm4_pcode,adm4_name,name_latin`. There is no `ref` or `fclass`; `source` and ADM4 fields are entirely empty.
- IDs are unique strings: 447,425 `way/...`, 2,617 `node/...`, 5 `relation/...`; maximum length 17. Preserve them as strings and add a stable part index when line normalization explodes a source geometry.
- Name coverage: 10,290 records have at least one of `name`, `name_en`, `name_am`, `name_latin`; 439,821 have no generic `name`. Keep source-name fields distinct from generated capital context.
- Surface is populated on 47,886 records (10.64%), so it is present but too incomplete to equate `is_main` with “all-weather” without a sourced rule.
- The approved main allowlist (`motorway`, `trunk`, `primary`, `secondary` and their `_link` classes) selects 14,581 line records. The approved “any road” allowlist adds `tertiary`, `tertiary_link`, `unclassified`, `residential`, `living_street`, `service`, `track` and `road`, selecting 395,963 lines. Construction, proposed, non-motor and ambiguous classes are excluded.
- Complete observed class counts (all geometries): `bridleway` 29, `bus_stop` 417, `construction` 1,869, `corridor` 1, `crossing` 1,761, `cycleway` 37, `emergency_access_point` 1, `emergency_bay` 1, `escape` 1, `footway` 4,974, `living_street` 271, `milestone` 2, `mini_roundabout` 29, `motorway` 128, `motorway_junction` 38, `motorway_link` 113, `path` 43,892, `pedestrian` 73, `platform` 1, `primary` 4,580, `primary_link` 646, `proposed` 1, `raceway` 2, `residential` 243,959, `rest_area` 37, `road` 147, `secondary` 5,449, `secondary_link` 309, `service` 17,632, `services` 21, `steps` 215, `stop` 1, `street_lamp` 57, `tertiary` 14,661, `tertiary_link` 235, `toll_gantry` 6, `track` 40,337, `traffic_mirror` 2, `traffic_signals` 84, `trunk` 3,165, `trunk_link` 355, `turning_circle` 165, `turning_loop` 4, `unclassified` 64,338, `yes` 1.

## 6. EthioSIS sample CSV

- 52,055 raw rows. Coordinates: 45 rows with blank/invalid pairs, 0 valid pairs outside the canonical extent, 52,004 unique valid coordinate pairs, 6 duplicate-coordinate groups (6 excess rows).
- Exact headers: `Plot Code`, `lon`, `lat`, `B (mg/kg)`, `Ca++ (mg/kg)`, `OC (%)`, `Co mg/kg`, `Cu mg/kg`, `Fe (mg/kg)`, `Mg (mg/kg)`, `Mn (mg/kg)`, `Molybdenium (mg/kg)`, `N total (%)`, `P (mg/kg)`, `K (mg/kg)`, `Silica (mg/kg)`, `Na (mg/kg)`, `S (mg/kg)`, `Zn (mg/kg)`, `Electrical Conductivity Property (dS/m)`, `pH`. Header-based mapping is mandatory; do not map columns by position.
- The raw CSV remains read-only. Its OC, N, P, pH, EC, boron and cobalt columns are not equivalent to bands present in the IDW file. In particular, raw `N total (%)` is not an application nitrogen source.

| Column | Populated | Blank | Nonpositive | Min | Median | Max | Positive p99.9 | Inventory note |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `B (mg/kg)` | 51,484 | 571 | 50 | -0.17 | 0.36 | 138.822 | 11.286 |  |
| `Ca++ (mg/kg)` | 52,032 | 23 | 0 | 24 | 3,458 | 37,891 | 26,131.814 |  |
| `OC (%)` | 6,700 | 45,355 | 0 | 0.325581 | 2.07 | 6.047 | 4.827 | not present in supplied IDW raster |
| `Co mg/kg` | 50,846 | 1,209 | 0 | 0.0005 | 1.996 | 18.71 | 14.606 |  |
| `Cu mg/kg` | 52,021 | 34 | 0 | 0.0025 | 2.73 | 195 | 13.299 |  |
| `Fe (mg/kg)` | 52,032 | 23 | 0 | 0.034811 | 121 | 822.166 | 565.9 |  |
| `Mg (mg/kg)` | 52,029 | 26 | 0 | 4 | 487 | 6,094 | 3,512.916 |  |
| `Mn (mg/kg)` | 52,031 | 24 | 0 | 1 | 143 | 704.465 | 543.985 |  |
| `Molybdenium (mg/kg)` | 44,114 | 7,941 | 0 | 0.000157 | 1.014 | 1,167.3 | 939.2 |  |
| `N total (%)` | 6,700 | 45,355 | 0 | 0.005 | 0.17 | 0.910153 | 0.522763 | retained raw only; excluded from application nitrogen |
| `P (mg/kg)` | 51,742 | 313 | 0 | 0.017 | 7 | 1,501 | 546.777 | not present in supplied IDW raster |
| `K (mg/kg)` | 52,031 | 24 | 0 | 9 | 365 | 4,368.73 | 2,489.61 |  |
| `Silica (mg/kg)` | 50,324 | 1,731 | 0 | 0.1 | 419 | 1,681 | 1,246.189 |  |
| `Na (mg/kg)` | 52,032 | 23 | 0 | 0.3 | 34.2 | 18,195.7 | 5,838.501 |  |
| `S (mg/kg)` | 51,949 | 106 | 0 | 0.25 | 11 | 10,772.9 | 7,425.646 |  |
| `Zn (mg/kg)` | 51,962 | 93 | 0 | 0.002739 | 1.8 | 375.242 | 41.204 |  |
| `Electrical Conductivity Property (dS/m)` | 37,101 | 14,954 | 0 | 0.017 | 90 | 16,560 | 2,779 | median 90 dS/m is suspicious; not used without unit/data review |
| `pH` | 38,602 | 13,453 | 0 | 2.26 | 6.39 | 65.1 | 9.432 | max 65.1 is implausible; clean per verified rule before loading |

The specified loader rule blanks nonpositive values and values above each column’s positive 99.9th percentile, then averages identical coordinates. The inventory quantifies the raw inputs but does not execute or approve a particular cleaning order. Loader reports must account for invalid coordinates, per-column blanks/outliers, and six coordinate merges.

## 7. Other input/configuration files

- `rainfall_seasonal_1km.json`: 9 records; band names, seasons, month sets and included years match the raster.
- Roads `metadata.json`, `README.txt`, `config.yaml`: agree on source, snapshot, licence and one roads export; metadata confirms mixed geometry and 450,047 rows.
- Seven `.aux.xml` files are PAM statistics sidecars only. They are not independent data sources and were not edited.
- Current scoring/configuration files outside `data/inputs`: `soil_thresholds.csv` (17 rows), `crops.csv` (6 rows), and misspelled `crop_requirments.csv` (36 rows). None contains literal `TBD`; `coarse_fragments` is explicitly `PROVISIONAL`, and several threshold notes still say units are unconfirmed. They belong to Phase 2 and must not reach a published map unreviewed.
- `data/scoring/soil_tv2.csv`, referenced in repository instructions as a prior mismatch, is not present in the current working tree. At audit time `soil_thresholds.csv` is the only threshold file on disk; canonical Phase 2 ownership still requires explicit confirmation before scoring.

## 8. Spot checks

| Point | DEM 1 km | DEM 90 m | Native slope 1 km | Native slope 90 m | Kiremt mean | Result |
|---|---:|---:|---:|---:|---:|---|
| Addis Ababa | 2,412.551 m | 2,410.776 m | 9.378% | 9.642% | 878.422 mm | plausible |
| Danakil/Dallol | -97.358 m | -93.633 m | 2.91% | 0.513456% | 94.099 mm | plausible |
| Bahir Dar | 1,792.666 m | 1,796.489 m | 1.83% | 0.69177% | 1,201.135 mm | plausible |
| Somali lowland | 513.683 m | 515.906 m | 0.393387% | 0.167125% | 23.51 mm | lower rainfall than highland as expected |

- Addis Ababa sampled about 2,413 m at 1 km and 2,411 m at 90 m, close to the rough 2,350 m sanity target given coordinate/product differences.
- Danakil sampled about -97 m at 1 km and -94 m at 90 m, confirming below-sea-level terrain.
- Bahir Dar kiremt mean is about 1,201 mm versus about 23.5 mm at the inspected Somali lowland point.
- Albers slope samples were transformed from EPSG:4326 before reading. For illustration only, Addis 90 m 9.642% converts to 5.507°; this exact sample conversion is not the common-grid mean-degree calculation.

## 9. Planned Phase 1 stack after sign-off

The core v1 plan has 38 bands: mask/admin 2, SoilGrids/Africa product 9, EthioSIS 12, terrain 2, land cover 1, rainfall 9 and temperature 3. The three optional access bands would make 41. No EthioSIS OC/N/P bands exist or will be synthesized.

| Input | Disposition |
|---|---|
| SoilGrids bands except phosphorus | copy after exact-grid checks; mask outside Ethiopia; use documented units |
| `sg_n` | copy in mg/kg without rescaling; sole nitrogen input |
| `sg_p` | approved AfricanSoilGrids input; exact-grid masked identity-copy in mg/100kg, with no scaling or back-transformation |
| EthioSIS IDW | approved: window-copy the first 1812 columns without interpolation; preserve nutrient-specific nulls independently of positive confidence |
| 90 m slope | derive percent-to-degree pixels under `data/interim/`, then average to canonical grid; preserve original |
| Existing 1 km slope | reference/validation only; do not relabel or convert its mean percentage as mean degrees |
| 1 km DEM, land cover, seasonal rainfall, temperature | consume after storage policy; retain/report inside-mask gaps and source provenance |
| Monthly rainfall | upstream provenance only; no rebuild |
| Admin polygons | approved as the full-resolution, publishable v04/2025-01-01 source; rasterize with deterministic ID mapping |
| Roads/capitals | approved motor-road classes recorded; PostGIS/PMTiles work deferred until database/tooling is supplied |
| Markets | unavailable pending actual market locations; not a core stack prerequisite |

## 10. Mismatches that must not be papered over

- Specifications requested mixed uint8/uint16/float32 bands and per-kind nodata in one GeoTIFF, while GeoTIFF storage has one physical datatype and generally one dataset nodata tag. `grid.yaml` says -9999; continuous-band docs say NaN. Mixed categorical/continuous overview resampling also needs a prototype.
- EthioSIS is one column wider than the grid and has nutrient-specific missingness despite positive confidence.
- Slope inputs are percentages, not degrees. Existing `slope_1km.tif` reaches 138.584% inside Ethiopia and cannot be inserted as `slope_mean_deg`.
- Phosphorus has no embedded unit/source metadata or band description. Its AfricanSoilGrids source, mg/100kg unit and no-transform policy are owner-confirmed, not inferred from its containing SoilGrids folder or approximate PAM statistics.
- Roads have no places layer, `ref` or `fclass`; mixed non-line features cannot be loaded as road lines, and the `surface` field is sparse rather than absent.
- Woreda polygons do not contain markets or populated centre fields. Gazetteer/administrative-capital points are not verified markets.
- Administrative geometries have version/date metadata but no sufficient licence/source/full-resolution declaration in the supplied files.
- The installed runtime/tooling is incomplete and no PostGIS target is configured.

## 11. Owner decisions and infrastructure deferral

| Gate | Owner decision (2026-11-23, with dated follow-ups) | Build consequence |
|---|---|---|
| Physical stack storage/nodata | Approved uniform float32 physical storage and canonical `-9999`; retain logical types in the manifest | Prototype and validate average versus mode overviews before assembling real bands; do not change `grid.yaml`. |
| EthioSIS width | Approved a window-copy of the first 1812 columns with no interpolation | Verify the discarded eastern column on every run and write an audit report. |
| EthioSIS confidence | Approved nutrient-specific nulls even where confidence is positive | Validate each nutrient independently; never impute missing nutrients. |
| Phosphorus | Final follow-up 2026-11-24: AfricanSoilGrids, not global SoilGrids; supplied mg/100kg values are approved as-is and require no scaling or change | **Scientific gate closed.** Include the lossless, country-masked `sg_p`; do not rescale, log/back-transform, clip or impute its valid values. Scoring remains Phase 2 work. |
| Admin boundaries | Approved the supplied v04/2025-01-01 polygons as full-resolution lookup/rasterization data with publication rights | Build the country mask and deterministic woreda-ID mapping from these files. |
| Road classes | Approved an explicit motor-road allowlist: main classes plus links; “any” also includes tertiary/link, unclassified, residential, living street, service, track and road | Use these classes when named-object work resumes; main does not mean all-weather. |
| Optional access bands | Deferred all three from core v1 | Add only through a later append-only stack version after road/centre validation. |
| Runtime | Approved a Python 3.12 project `.venv`, reusing system GDAL and installing pinned missing packages only inside it | Add the dependency manifest and run implementation/tests through `.venv`. |
| PostGIS/tooling | Deferred named-object and PMTiles work because no test database or toolchain is available | Implement the raster foundation first; do not create untestable loaders/assets now. |
| Markets | Actual locations remain pending | Keep unavailable; no centroid/capital proxy. |

Phase 2-only questions (threshold canonicality, provisional scientific values and crop/scoring review) do not authorize Phase 1 scoring work and should not delay non-scoring inventory sign-off.

## 12. Sign-off record

The decisions in section 11 were explicitly recorded through the implementation conversation. They authorize the listed raster-foundation work, including final stack assembly after the owner's AfricanSoilGrids/as-is confirmation in section 14. No source data was changed. Named-object/PostGIS/PMTiles infrastructure is still deferred; this is not complete Phase 1 sign-off.

## 13. Initial phosphorus unit confirmation and interim staging (2026-11-24)

This records the earlier unit-only stage. The subsequent owner approval in section 14 supersedes its pending publication status and answers the questions below.

### Confirmed and verified

- The owner identified the phosphorus input under `data/inputs/SoilGrids/` and confirmed **mg/100kg**. The inspected file is `data/inputs/SoilGrids/p_0-30cm_mean_1km_4326.tif`.
- It is a single Float32 band, EPSG:4326, pixel-is-area, with canonical shape/affine and `-9999` nodata. `config.grid.assert_matches_grid` passes without a warp or crop.
- Band description and unit tags are empty; embedded scale and offset are absent. The explicit singleton filename alias supplies `sg_p`; it is not an assumed multiband index or proof of upstream provenance. The PAM sidecar contains approximate statistics only, so reported values below are fresh full scans.
- Source SHA-256: `0d2e8c361bf6221b9255ca3dc0f5f7a077aeb67e8d5e74809f4a99e7a2557ab0`.

| Inside Ethiopia | Result |
|---|---:|
| Country-mask cells | 1,338,026 |
| Valid phosphorus cells | 1,332,610 |
| Preserved missing cells | 5,416 |
| Minimum (mg/100kg) | 174.62637329101562 |
| Median (mg/100kg) | 1,596.827392578125 |
| Maximum (mg/100kg) | 9,091.994140625 |
| Cells below the provisional specification limit of 70 mg/100kg | 0 |
| Cells above the provisional specification limit of 4,000 mg/100kg | 2,200 |

The provisional specification limits are diagnostics only, not accepted scientific bounds or scoring thresholds. No values were clipped. For dimensional reference only, dividing mg/100kg by 100 gives mg/kg: the minimum/median/maximum are equivalent to approximately 1.746/15.968/90.920 mg/kg. **The raster remains mg/100kg and was not rescaled.**

### Local artifacts and checks

- Script: `data/interim/phosphorus/stage_phosphorus.py`
- Interim-only raster: `data/interim/phosphorus/phosphorus_staged_pending.tif`
- Audit: `data/interim/phosphorus/phosphorus_report.json`
- Tests: `data/interim/phosphorus/test_stage_phosphorus.py` (20 synthetic checks passing).
- Every valid inside-mask value round-trips bit-for-bit; every outside-mask cell is canonical nodata. Source gaps are preserved, with no resampling, extra scaling, log/back-transform, clipping or imputation.
- Raster and report are byte-identical on rerun. SHA-256 checks for the source, its PAM sidecar and the mask are unchanged; all 219 input files retain their sizes/modification times and the file inventory is unchanged.
- The raster and report explicitly say `stack_eligible: false` and `publication_ready: false`. Nothing was copied to `data/processed/`; this is not a final COG/manifest delivery.

### Questions at the initial unit-only stage (resolved in section 14)

1. What is the exact upstream phosphorus product/version/source link and measurement meaning (for example, the extractable-P method)? Being in the SoilGrids folder does not establish that it is an ISRIC SoilGrids product; the specification identifies a separate Africa-wide product.
2. Are the supplied pixel values already physical mg/100kg after any product-specific scaling/back-transformation? Confirm the upstream processing, or provide the required conversion if they remain encoded. No processing history or nonlinear transform was inferred from the observed range.

These questions initially held `P1-02`/`P1-11` open. The later AfricanSoilGrids/as-is confirmation closes that gate; PostGIS loaders and PMTiles remain separately deferred for lack of the approved test database/toolchain.

## 14. Final AfricanSoilGrids phosphorus approval (2026-11-24)

- **Owner-confirmed product:** AfricanSoilGrids, not the normal/global SoilGrids product.
- **Owner-confirmed unit:** mg/100kg.
- **Owner-approved operation:** use the supplied pixel values as-is; no scaling, unit conversion, log transform or back-transformation is required or applied. The standard Ethiopia mask and source nodata handling remain part of the approved stack policy.
- **Approval scope:** the existing `data/inputs/SoilGrids/p_0-30cm_mean_1km_4326.tif`, whose SHA-256 remains the value recorded in section 13. Product version/source URL and upstream processing history are not invented; the owner's explicit as-is approval is the recorded provenance authority.
- **Approved artifact:** `data/interim/phosphorus/phosphorus_staged.tif`; the earlier `phosphorus_staged_pending.tif` is a superseded local artifact and is not consumed by the stack build.
- **Current audit:** `data/interim/phosphorus/phosphorus_report.json`, now `status: passed`, `stack_eligible: true`, with the AfricanSoilGrids identity and owner approval recorded. No publication blocker remains for this input.
- **Checks:** 20 regression tests pass; every valid phosphorus pixel round-trips bit-for-bit; the 5,416 inside-country gaps are preserved; source/mask hashes are unchanged; reruns produce byte-identical approved raster and report.
- **Range findings:** the 2,200 cells above the specification's provisional diagnostic limit are retained and reported, not treated as scientific errors or scored/clipped. No provisional scoring threshold is introduced.

`P1-02` and `P1-11` are complete. All 38 core-band ingredients exist locally and final raster assembly/acceptance can proceed. This does not authorize installing named-object infrastructure or starting scoring/API/UI work.
