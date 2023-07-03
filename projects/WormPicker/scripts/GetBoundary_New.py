import csv
import subprocess
import sys
from optparse import OptionParser, OptionGroup
import gzip
import os
import argparse
import ast
import math
from module_crs_converter import trans4mEPSG
import re

#########################################################   HELP   #########################################################################
usage="python %prog --source input.csv "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--LayerInfoFile", dest="INFO", help="The Layer Info file.")


parser.add_option_group(group)
(options, args) = parser.parse_args()


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


layer_names=[]
min_x_values = []
max_x_values = []
min_y_values = []
max_y_values = []


#print(layer_names)
with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/new_layer_info_WCS.csv", newline='') as csvfile:
#with open(options.INFO, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    #next(reader) 
    for row in reader:
        layer=row["CoverageID"]
        crs=row['CRS']
        min_x=float(row['minlat'])
        max_x=(float(row['maxlat']))
        min_y=(float(row['minlong']))
        max_y=(float(row['maxlong']))
        match= re.search("EPSG",crs)
        if crs=="EPSG/0/4326":
            min_y,min_x=trans4mEPSG("EPSG:4326","EPSG:3035", min_y, min_x)
            max_y,max_x=trans4mEPSG("EPSG:4326","EPSG:3035",max_y, max_x)
            min_x_values.append(min_x)
            max_x_values.append(max_x)
            min_y_values.append(min_y)
            max_y_values.append(max_y)
            layer_names.append(layer)
        elif match:
            crs_indic = crs.replace("EPSG/0/", "EPSG:")
            min_y, min_x = trans4mEPSG(crs_indic,"EPSG:3035", min_y, min_x)
            max_y, max_x = trans4mEPSG(crs_indic,"EPSG:3035",max_y, max_x)
            min_x_values.append(float(row['minlat']))
            max_x_values.append(float(row['maxlat']))
            min_y_values.append(float(row['minlong']))
            max_y_values.append(float(row['maxlong']))
            layer_names.append(layer)
        else:
            next
    #print(layer_names)
    #print(min_x_values)

# Calculate the common boundary box
common_min_x = min([x for x in min_x_values if x != float('inf')])
common_max_x = max([x for x in max_x_values if x != float('inf')])
common_min_y = min([x for x in min_y_values if x != float('inf')])
common_max_y = max([x for x in max_y_values if x != float('inf')])

# Print the common boundary box
print(f"{common_min_x}:{common_max_x}:{common_min_y}:{common_max_y}")
