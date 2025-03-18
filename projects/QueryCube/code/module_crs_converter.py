import math
from osgeo import gdal, osr

def trans4mEPSG(InputCRS,OutputCRS,y,x):
    src_crs = InputCRS
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(src_crs)
    tgt_crs = OutputCRS
    tgt_srs = osr.SpatialReference()
    tgt_srs.SetFromUserInput(tgt_crs)
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    x_t, y_t, z_t = transform.TransformPoint(x, y, 0)
    return x_t, y_t

trans4mEPSG("EPSG:3035","EPSG:4258",-64237.727615634445, 8393681.506823568)
trans4mEPSG("EPSG:4258","EPSG:3035",7.05332, 46.93947)
trans4mEPSG("EPSG:4258","EPSG:3035", 46.93947,7.05332)
trans4mEPSG("EPSG:4326","EPSG:3035",46.8136889,13.50794792)
trans4mEPSG("EPSG:3035","EPSG:4326",2651789.8737778864, 4096504.5378351957)
