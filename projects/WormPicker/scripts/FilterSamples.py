import sys
import csv
from optparse import OptionParser, OptionGroup
import gzip


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

#lat = float(options.lat)
lat_l = float(str.split(options.bound,":")[0])
lat_u = float(str.split(options.bound,":")[1])
lon_l =float(str.split(options.bound,":")[2])
lon_u = float(str.split(options.bound,":")[3])


#lon_l = float(str.split(options.bound, ":")[0]) 
#lon_u = float(str.split(options.bound, ":")[1]) 

with open(options.source, 'r') as f:
    reader = csv.reader(f)
    header=next(reader) ##added
    data = list(reader)

filtered_data=[]

for row in data:
    if row[3] == "NA":
        continue
    elif float(row[3]) > lat_l and float(row[4]) < lat_u and float(row[4]) > lon_l and float(row[4]) < lon_u:
        filtered_data.append(row)

with open('output/coveredsamples_wcs.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header) ###added
    writer.writerows(filtered_data)