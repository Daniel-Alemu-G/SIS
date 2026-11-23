from pathlib import Path
import yaml
from affine import Affine

_cfg = yaml.safe_load(
    (Path(__file__).parents[1] / "pipeline" / "grid.yaml").read_text()
)

CRS = _cfg["crs"]
CPD = _cfg["cells_per_degree"]
RES = 1 / CPD
WEST, NORTH = _cfg["origin"]["lon"], _cfg["origin"]["lat"]
COLS, ROWS = _cfg["size"]["cols"], _cfg["size"]["rows"]
NODATA = _cfg["nodata"]

TRANSFORM = Affine(RES, 0, WEST, 0, -RES, NORTH)
EAST = WEST + COLS * RES    # 48.0
SOUTH = NORTH - ROWS * RES  # 3.4

PROFILE = dict(
    driver="GTiff", crs=CRS, transform=TRANSFORM,
    width=COLS, height=ROWS, nodata=NODATA,
    compress="deflate", tiled=True,
)