import math
from os import PathLike
from pathlib import Path

import yaml
from affine import Affine

_cfg = yaml.safe_load(
    (Path(__file__).parents[1] / "pipeline" / "grid.yaml").read_text()
)

CRS = _cfg["crs"]
CPD = _cfg["cells_per_degree"]
RES = float(_cfg["cell_size_deg"])
WEST, NORTH = _cfg["origin"]["lon"], _cfg["origin"]["lat"]
COLS, ROWS = _cfg["size"]["cols"], _cfg["size"]["rows"]
NODATA = _cfg["nodata"]

TRANSFORM = Affine(*_cfg["transform"])
EAST = _cfg["extent"]["east"]
SOUTH = _cfg["extent"]["south"]

PROFILE = dict(
    driver="GTiff", crs=CRS, transform=TRANSFORM,
    width=COLS, height=ROWS, nodata=NODATA,
    compress="deflate", tiled=True,
)


class GridMismatchError(ValueError):
    """Raised when a common-grid raster does not match ``grid.yaml``."""


def _validate_grid_config() -> None:
    """Fail at import if the canonical YAML contradicts itself."""
    checks = {
        "cell_size_deg": (RES, 1 / CPD),
        "transform.a": (TRANSFORM.a, RES),
        "transform.b": (TRANSFORM.b, 0),
        "transform.c": (TRANSFORM.c, WEST),
        "transform.d": (TRANSFORM.d, 0),
        "transform.e": (TRANSFORM.e, -RES),
        "transform.f": (TRANSFORM.f, NORTH),
        "extent.east": (EAST, WEST + COLS * RES),
        "extent.south": (SOUTH, NORTH - ROWS * RES),
    }
    inconsistent = [
        f"{name}: {actual!r} != {expected!r}"
        for name, (actual, expected) in checks.items()
        if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12)
    ]
    if inconsistent:
        raise ValueError(
            "pipeline/grid.yaml is internally inconsistent: "
            + "; ".join(inconsistent)
        )


def grid_mismatches(dataset_or_path) -> list[str]:
    """Return common-grid mismatches for a GDAL raster dataset or path.

    The comparison checks equivalent CRS definitions, dimensions, every
    geotransform component (including rotation), and pixel-is-area semantics.
    It never repairs, crops, warps, or opens a path in update mode.
    """
    from osgeo import gdal, osr

    opened_here = isinstance(dataset_or_path, (str, PathLike))
    if opened_here:
        dataset = gdal.OpenEx(
            str(dataset_or_path), gdal.OF_RASTER | gdal.OF_READONLY
        )
        if dataset is None:
            return [f"cannot open raster read-only: {dataset_or_path}"]
        label = str(dataset_or_path)
    elif hasattr(dataset_or_path, "GetGeoTransform"):
        dataset = dataset_or_path
        label = dataset.GetDescription() or "<GDAL dataset>"
    else:
        raise TypeError("expected a raster path or GDAL Dataset")

    mismatches = []
    if dataset.RasterXSize != COLS or dataset.RasterYSize != ROWS:
        mismatches.append(
            "shape: "
            f"{dataset.RasterXSize}x{dataset.RasterYSize} != {COLS}x{ROWS}"
        )

    actual_transform = dataset.GetGeoTransform(can_return_null=True)
    expected_transform = TRANSFORM.to_gdal()
    if actual_transform is None:
        mismatches.append("transform: missing")
    else:
        for index, (actual, expected) in enumerate(
            zip(actual_transform, expected_transform, strict=True)
        ):
            if not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
                mismatches.append(
                    f"transform[{index}]: {actual!r} != {expected!r}"
                )

    actual_crs = dataset.GetSpatialRef()
    if actual_crs is None:
        mismatches.append("crs: missing")
    else:
        expected_crs = osr.SpatialReference()
        expected_crs.SetFromUserInput(CRS)
        actual_crs = actual_crs.Clone()
        actual_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        expected_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        if not actual_crs.IsSame(expected_crs):
            mismatches.append(
                f"crs: {actual_crs.GetName()!r} is not equivalent to {CRS}"
            )

    area_or_point = dataset.GetMetadataItem("AREA_OR_POINT")
    if (area_or_point or "").lower() != "area":
        mismatches.append(
            f"pixel convention: AREA_OR_POINT={area_or_point!r}, expected 'Area'"
        )

    if opened_here:
        dataset = None
    return [f"{label}: {item}" for item in mismatches]


def assert_matches_grid(dataset_or_path) -> None:
    """Raise :class:`GridMismatchError` unless a raster matches the grid."""
    mismatches = grid_mismatches(dataset_or_path)
    if mismatches:
        raise GridMismatchError("common-grid mismatch:\n- " + "\n- ".join(mismatches))


_validate_grid_config()
