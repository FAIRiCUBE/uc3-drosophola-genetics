import sys
import csv
from optparse import OptionParser, OptionGroup
import gzip
import math

#########################################################   HELP   #########################################################################
usage="python %prog --source input.csv --boundary-lat lat --bundary-lon lon > covered.csv"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--source", dest="source", help="The source CSV with lat and lon in the 4th and 5th column.")
#parser.add_option("--boundary-lat", dest="lat", help="")
parser.add_option("--boundary", dest="bound", help="")


parser.add_option_group(group)
(options, args) = parser.parse_args()


def convert_4326_to_3035(longitude, latitude):
    # Convert degrees to radians
    lon_rad = math.radians(longitude)
    lat_rad = math.radians(latitude)
    # ETRS89 (ETRS-LAEA) projection parameters
    central_meridian = math.radians(10.0)
    scale_factor = 1.0
    false_easting = 4321000.0
    false_northing = 3210000.0
    # Lambert Azimuthal Equal Area projection formulas
    rho = 6378137.0 * scale_factor
    x = rho * math.cos(lat_rad) * math.sin(lon_rad - central_meridian) + false_easting
    y = (rho * (math.cos(lat_rad) * math.cos(lon_rad - central_meridian))-(rho * math.sin(lat_rad)) + false_northing)
    return x, y

#lat = float(options.lat)
lat_l = float(str.split(options.bound,":")[0])
lat_u = float(str.split(options.bound,":")[1])
lon_l =float(str.split(options.bound,":")[2])
lon_u = float(str.split(options.bound,":")[3])


#lon_l = float(str.split(options.bound, ":")[0]) 
#lon_u = float(str.split(options.bound, ":")[1])
lat_l=11541132.485614082
lat_u=-3168137.0
lon_l=6763968.203153217
lon_u=4321000.0

filtered_data=[]

with open(options.source, 'r') as f:
    reader = csv.DictReader(f)
    header = next(reader)
    for row in reader:
        lat=row["lat"]
        long=row["long"]
        if lat and long == "NA":
            continue
        else:
            long,lat=convert_4326_to_3035(float(long),float(lat))
            #print(long,lat)
            if float(lat) < lat_l and float(lat) > lat_u and float(long) < lon_l and float(long) > lon_u:
                #print("YES")
                filtered_data.append(row.values())

print(filtered_data)

with open('output/coveredsamples_wcs.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(filtered_data)