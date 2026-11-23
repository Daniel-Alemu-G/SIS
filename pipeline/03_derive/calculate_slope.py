from osgeo import gdal

gdal.UseExceptions()

def process_accurate_slope(input_dem, output_metric_dem, output_slope, res=90):
    print("Warping to metric projection...")
    gdal.Warp(
        output_metric_dem, input_dem,
        dstSRS='+proj=aea +lat_1=3 +lat_2=15 +lat_0=9 +lon_0=40 +datum=WGS84 +units=m +no_defs',
        resampleAlg='bilinear',
        xRes=res, yRes=res,
        srcNodata=-32768, dstNodata=-32768,
        creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']
    )

    print("Calculating slope...")
    gdal.DEMProcessing(
        output_slope, output_metric_dem, 'slope',
        format='GTiff',
        computeEdges=True,
        slopeFormat='percent',
        alg='Horn',
        creationOptions=['COMPRESS=DEFLATE', 'TILED=YES']
    )
    print(f"[ok] Slope saved to: {output_slope}")

process_accurate_slope('dem_glo90.tif', 'dem_albers_90m.tif', 'slope_albers_90m.tif')