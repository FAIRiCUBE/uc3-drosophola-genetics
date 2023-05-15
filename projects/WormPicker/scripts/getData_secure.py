import csv
import subprocess
import sys
from optparse import OptionParser, OptionGroup
import gzip
import os


#########################################################   HELP   #########################################################################
usage="python %prog --source input.csv "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--source", dest="source", help="The source CSV with the names of the sample information.")

parser.add_option_group(group)
(options, args) = parser.parse_args()



url = 'https://fairicube.rasdaman.com/rasdaman/ows'
headers = {"Content-Type": "application/x-www-form-urlencoded"}

firstline="sampleId,dlt,esm,ft,imp,tcd,ww"


with open(options.source, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader) 
    for row in reader:
        lat = row[3]
        lon = row[4]
        if lat != 'NA' and lon != 'NA':
            lat2 = round((float(lat)+0.0001),5)
            lon2 = round((float(lon)+0.0001),5)
            lat_range = f"Lat({lat}%3A{lat2})"
            lon_range = f"Lon({lon}%3A{lon2})"
            query = "for%20%24dlt%20in%20(dominant_leaf_type)%2C%20%24esm%20in%20(european_settlement_map)%2C%20%24ft%20in%20(high_resolution_layer_forest_type)%2C%20%24imp%20in%20(imperviousness)%2C%20%24tcd%20in%20(tree_cover_density)%2C%20%24ww%20in%20(water_and_wetness)%20return%20encode((%20%7Bimp%3A%24imp%3B%20dlt%3A%24dlt%5Bansi(%222018-01-01T08%3A12%3A02.000Z%22)%5D%3B%20esm%3A%20%24esm%5B%20ansi(%222018-01-01T08%3A12%3A02.000Z%22)%20%5D%3B%20ft%3A%24ft%5B%20ansi(%222018-01-01T08%3A12%3A02.000Z%22)%5D%3B%20tcd%3A%24tcd%5B%20ansi(%222018-01-01T08%3A12%3A02.000Z%22)%5D%3B%20ww%3A%24ww%5B%20ansi(%222018-01-01T08%3A12%3A02.000Z%22)%5D%7D)%20%5B"+f"Lat({lat}%3A{lat2})%2C%20Lon({lon}%3A{lon2})%5D%20%2C%20%22text%2Fcsv%22)"
            data = "SERVICE=WCS&VERSION=2.0.1&REQUEST=ProcessCoverages&QUERY="
            samplename=str({row[0]})
            output_file='/media/ssteindl/fairicube/rasda/output/MetaForCoveredSamples.csv'
            result = subprocess.run(['curl',
                    '-u', 'XXXXXXXXXXXXXXXX:XXXXXXX',
                    '-X', 'POST',
                    '-H', 'Content-Type: ' + headers["Content-Type"],
                    '-d', '"' + data + query + '"',
                    url],
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
            res1= result.stdout.decode('utf-8').strip()[2:-2]


            with open(output_file, 'a') as f:
                    f.write(firstline) ##added April 3rd 
                    f.write(str(samplename)[2:-2] + ',' + res1.replace(" ", ",") + '\n')
        else:
            print("NANA")




