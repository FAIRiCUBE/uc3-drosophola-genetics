import sys
import csv
from optparse import OptionParser, OptionGroup
import gzip
import math
import os
import module_crs_converter

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

print(os.getcwd())


#lat = float(options.lat)
lat_l = float(str.split(options.bound,":")[0])          #latitude lower boundary
lat_u = float(str.split(options.bound,":")[1])          #latitute upper boundary
lon_l =float(str.split(options.bound,":")[2])           #longitude lower boundary
lon_u = float(str.split(options.bound,":")[3])          #longitude upper boundary


#lat_l = float(str.split(output,":")[0])          #latitude lower boundary
#lat_u = float(str.split(output,":")[1])          #latitute upper boundary
#lon_l =float(str.split(output,":")[2])           #longitude lower boundary
#lon_u = float(str.split(output,":")[3])          #longitude upper boundary



#lon_l = float(str.split(options.bound, ":")[0]) 
#lon_u = float(str.split(options.bound, ":")[1])
#lat_l=11541132.485614082
#lat_u=-3168137.0
#lon_l=6763968.203153217
#lon_u=4321000.0

filtered_data=[]

with open(options.source, 'r') as f:
#with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv", 'r') as f:
    reader = csv.DictReader(f)
    header = next(reader)
    for row in reader:
            lat=row["lat"]
            long=row["long"]
            if lat == 'NA' and long == 'NA':
                continue
            else:
                print(long,lat)
            long,lat=module_crs_converter.trans4mEPSG("EPSG:4326","EPSG:3035",float(long),float(lat))
            #print(long,lat)
            if float(lat) > lat_l and float(lat) < lat_u and float(long) > lon_l and float(long) < lon_u:
                #print("YES")
                filtered_data.append(row.values())

print(filtered_data)

with open('output/coveredsamples_wcs.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(filtered_data)