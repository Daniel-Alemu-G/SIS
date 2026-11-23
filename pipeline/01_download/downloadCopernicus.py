import os
import requests
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from pystac_client import Client

TILE_DIR = "dem_tiles"
OUT = "ethiopia_copernicus_dem_90m.tif"
os.makedirs(TILE_DIR, exist_ok=True)

# 1. Boundary
gdf = gpd.read_file("ethiopia.geojson").to_crs("EPSG:4326")
geom = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union

# 2. Find tiles
catalog = Client.open("https://earth-search.aws.element84.com/v1")
items = catalog.search(collections=["cop-dem-glo-90"], bbox=tuple(gdf.total_bounds)).items()

from shapely.geometry import shape
tiles = [i for i in items if shape(i.geometry).intersects(geom)]
print(f"{len(tiles)} tiles to download")

# 3. Download over HTTPS


paths = []
for item in tiles:
    href = item.assets["data"].href
    if href.startswith("s3://"):
        bucket, key = href[5:].split("/", 1)
        href = f"https://{bucket}.s3.amazonaws.com/{key}"
    path = os.path.join(TILE_DIR, os.path.basename(href))
    paths.append(path)
    if os.path.exists(path):
        continue
    print("Downloading", os.path.basename(href))
    with requests.get(href, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(path + ".part", "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    os.rename(path + ".part", path)

# 4. Merge to disk (windowed, low memory)
print("Merging...")
merge(
    paths,
    dst_path=OUT,
    dst_kwds={"compress": "deflate", "predictor": 3, "tiled": True, "BIGTIFF": "YES"},
)
print("Done:", OUT)