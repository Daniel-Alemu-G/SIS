#!/usr/bin/env python3
"""
Reproject and resample a Copernicus DEM GLO-90 raster to EPSG:4326 and resample to 1 km
(30 arc-seconds = 1/120 degree), snapped to the same global grid as the soil grids.

Usage:
    python resample_dem.py [IN_FILE] [-o OUT_FILE]
"""
import argparse
from pathlib import Path
import numpy as np
import rasterio
from rasterio.warp import reproject, transform_bounds, Resampling
from rasterio.transform import from_origin

RES = 1 / 120  # ~0.00833 deg, ~1 km at the equator
DST_CRS = "EPSG:4326"
NODATA = -9999.0


def process(src_path: Path, out_path: Path, method: Resampling):
    with rasterio.open(src_path) as src:
        left, bottom, right, top = transform_bounds(
            src.crs, DST_CRS, *src.bounds, densify_pts=101
        )
        # Snap outward to the 1/120 degree grid to align perfectly with other datasets
        left = np.floor(left / RES) * RES
        right = np.ceil(right / RES) * RES
        bottom = np.floor(bottom / RES) * RES
        top = np.ceil(top / RES) * RES

        width = int(round((right - left) / RES))
        height = int(round((top - bottom) / RES))
        transform = from_origin(left, top, RES, RES)

        dst = np.full((height, width), NODATA, dtype="float32")
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=transform,
            dst_crs=DST_CRS,
            dst_nodata=NODATA,
            resampling=method,
        )

        profile = src.profile.copy()
        profile.update(
            crs=DST_CRS, transform=transform, width=width, height=height,
            dtype="float32", nodata=NODATA, count=1,
            compress="deflate", predictor=3, tiled=True,
        )
        # Drop stale tiling keys if source was striped/tiled differently
        profile.pop("blockxsize", None)
        profile.pop("blockysize", None)
        profile["blockxsize"] = profile["blockysize"] = 256

    with rasterio.open(out_path, "w", **profile) as out:
        out.write(dst, 1)
    print(f"[ok] {src_path.name} -> {out_path.name} ({width}x{height})")


def resample():
    p = argparse.ArgumentParser()
    p.add_argument("in_file", type=Path, help="Path to the single Copernicus DEM GLO-90 TIFF file")
    p.add_argument("-o", "--out", type=Path, default=None, help="Optional output TIFF file path")
    p.add_argument("--method", default="bilinear",
                   choices=["average", "bilinear", "cubic", "med"],
                   help="resampling method (default: bilinear for elevation)")
    a = p.parse_args()

    if a.out:
        out_file = a.out
    else:
        out_file = a.in_file.parent / f"{a.in_file.stem}_1km_4326.tif"
        
    out_file.parent.mkdir(parents=True, exist_ok=True)
    method = getattr(Resampling, a.method)

    process(a.in_file, out_file, method)

if __name__ == "__main__":
    resample()