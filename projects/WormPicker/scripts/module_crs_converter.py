import math
from osgeo import gdal, osr

#def convert_4326_to_3035(longitude, latitude):
#    # Convert degrees to radians
#    lon_rad = math.radians(longitude)
#    lat_rad = math.radians(latitude)
#    # ETRS89 (ETRS-LAEA) projection parameters
#    central_meridian = math.radians(10.0)
#    scale_factor = 1.0
#    false_easting = 4321000.0
#    false_northing = 3210000.0
#    # Lambert Azimuthal Equal Area projection formulas
#    rho = 6378137.0 * scale_factor
#    x = rho * math.cos(lat_rad) * math.sin(lon_rad - central_meridian) + false_easting
#    y = (rho * (math.cos(lat_rad) * math.cos(lon_rad - central_meridian))-(rho * math.sin(lat_rad)) + false_northing)
#    return x, y


def trans4m2wgs84(InputCRS, x,y):
    src_crs = InputCRS
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(src_crs)
    tgt_crs = 'EPSG:4326' 
    tgt_srs = osr.SpatialReference()
    tgt_srs.SetFromUserInput(tgt_crs)
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    x_t, y_t, z_t = transform.TransformPoint(x, y, 0)
    return x_t, y_t


#trans4mwsg84('EPSG:3035',2639692.4010280264,4588860.382936129)

def trans4mfromwgs84(InputCRS, x,y):
    src_crs = 'EPSG:4326'
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(src_crs)
    tgt_crs = InputCRS
    tgt_srs = osr.SpatialReference()
    tgt_srs.SetFromUserInput(tgt_crs)
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    x_t, y_t, z_t = transform.TransformPoint(x, y, 0)
    return x_t, y_t

#trans4mfromwgs84('EPSG:3035',46.8136889,13.50794792)
