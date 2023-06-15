###


### 0. setUp and credentials

import requests
from dotenv import dotenv_values
import json
import xmltodict
import csv
from dateutil import parser     
from datetime import datetime
import math
from osgeo import gdal, osr

from module_crs_converter import trans4m2wgs84, trans4mfromwgs84

#####
#OPTIONS
usage="python %prog --source input.csv "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--source", dest="source", help="The source CSV with the names of the sample information.")
parser.add_option("--layers", dest="layers", help="The list of layers to get data from.")

parser.add_option_group(group)
(options, args) = parser.parse_args()

#####
env_vars = dotenv_values()

# Access specific environment variables
rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"

### 1. get a layer




def getlayer():
    with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/layer_info_WCS.csv", newline='') as csvfile:
        sample_reader = csv.reader(csvfile)
        row0 = next(sample_reader)
        row1 = next(sample_reader)
        row2 = next(sample_reader)
    return row2



trans4mEPSG("EPSG:4326","EPSG:3035",46.8136889,13.50794792)
#trans4mEPSG("EPSG:4326","EPSG:3035",13.50794792,46.8136889)
coords=trans4mEPSG("EPSG:4326","EPSG:3035",46.8136889,13.50794792)

latitude_request = coords[1]
longitude_request = coords[0]
resolution=int(layer[8])
crs_orig=layer[1]
crs = crs_orig.replace("EPSG/0/", "EPSG:")

if crs != "EPSG:4326":
    strlat="&subset=Y({},{})"
    strlon = "&subset=X({},{})"
else:
    strlat="&subset=Lat({},{})"
    strlon = "&subset=Long({},{})"

layer=getlayer()
coverage_id = layer[0]
Qs=((1,1),(-1,1),(-1,-1),(1,-1))

facts=[]
values=[]
resolutions=[]
cellnumbers=[]
for k in range (0,40):
    for i in range(0,4):
        x1= float(longitude_request)
        x2= x1 + resolution*Qs[i][0]*k*0.5*0.5
        y1= float(latitude_request)
        y2= y1 + resolution*Qs[i][1]*k*0.5*0.5
        if x1 > x2:
            z1=x1
            x1=x2
            x2=z1
        else:
            h=3
        if y1 > y2:
            z2=y1
            y1=y2
            y2=z2
        else:
            h=5
        ###fill in 
        subset_lat= strlon.format(y1,y2)
        subset_long = strlat.format(x1,x2)
        print(subset_lat, subset_long)
        #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
        request_cov_id = "&COVERAGEID=" + layer[0]
        request_encode_format = "&FORMAT=text/csv"
        #print(layer,lat,lon, row["minlat"], row["maxlat"])
        response = requests.get(base_wcs_url + "&request=GetCoverage" + request_cov_id + subset_lat + subset_long + request_encode_format,auth=(rasdaman_username, rasdaman_password), verify=False)
        value=response.text
        xi=value.split("},{")
        phi=(xi[0].split("', '"))
        cells=len(phi[0].split(','))
        facts.append(("Quadrant:",i,"Resolution", k*0.5*0.5, "Exclam:", len(xi), "Within:", cells))
        resolutions.append(k*0.5*0.5)
        cellnumbers.append(cells)
        values.append((k*0.5*0.5, i, value))
        
        #print(value)

### 2. set a starting point 

### 3. define directions that will be testet

### 4. test all and get results

### 5. table results 

### 6. find pattern in result 


