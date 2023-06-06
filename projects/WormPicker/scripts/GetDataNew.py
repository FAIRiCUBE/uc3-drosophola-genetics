import requests
from dotenv import dotenv_values
import json
import xmltodict
import csv
from dateutil import parser     
from datetime import datetime
import math
from osgeo import gdal, osr

import module_crs_converter

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

my_dict={} 
with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/coveredsamples_wcs.csv", newline='') as csvfile:
    reader_outer = csv.DictReader(csvfile)
    for s_row in reader_outer:
        latitude_original = s_row["lat"]
        longitude_original = s_row["long"]
        name_date = s_row["sampleId"]
        date = str(name_date[-10:])
        time_asked= date + 'T00:00:00Z'
        print("ORIGINAL COORDINATES OF SAMPLE\n:", name_date, latitude_original, longitude_original, "DATE:", time_asked)
        #datetime3 = parser.isoparse(time_asked)
        sample_result=[]
        #print(datetime3)
        layer=None
        converted_date=None
        my_dict[name_date]=[]  
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            wanted=parser.isoparse(converted_date)
            #print(wanted)
        except ValueError:
            #converted_date = date
            continue
        #print(converted_date)
        #print(lat)
        #print(date)
        #print(converted_date)
        with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/layer_info_WCS.csv", newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            #next(reader) 
            for row in reader:
                if any(value == "NA" for value in row.values()):
                    continue
                layer=row["CoverageID"]
                #print(layer)
                crs=row["CRS"]
                resolution=0
                lat = latitude_original
                lon = longitude_original
                time_min= row["f"]
                time_max= row["t"]
                print("LIMITS: ", time_min, time_max, "MI\n", row["minlat"], row["maxlat"] )
                datetime1 = parser.isoparse(time_min)
                datetime2 = parser.isoparse(time_max)
                #print(layer, crs, lat, lon, resolution)
                #print(converted_date, datetime1, datetime2) 
                #print(lat,lon)
                strlat="&subset=Y({},{})"
                strlon = "&subset=X({},{})"
                #print("SUBSETS DEFAULT",strlon,strlat)
                #print(abs(float(row["resolution_dim1"])))
                #print(datetime2, wanted, datetime1)
                if wanted > datetime1 and wanted < datetime2:
                    print("YES")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                else:
                    print("NOT COVERED TEMPORALLY")
                if abs(float(row["resolution_dim1"])) == abs(float(row["resolution_dim2"])):
                    resolution=abs(float(row["resolution_dim1"]))
                    #print("AT least axis resolution is identical")
                else:
                    #print("NOT same resolition in axis")
                    resolution=abs(float(row["resolution_dim2"])) ###still not sure why and what to do but spaceholder
                #print(resolution)
                #print(crs)
                if crs=="EPSG/0/4326":
                    lat=lat
                    lon=lon
                    strlat="&subset=Lat({},{})"
                    strlon = "&subset=Long({},{})"
                    #print(crs,lon,lat)
                elif crs=="EPSG/0/3035":
                    lon,lat= trans4mfromwgs84('EPSG:3035',float(lon), float(lat))
                    #print(crs,lon,lat)
                #elif crs=="EPSG/0/32631":
                    #lon,lat= convert_4326_to_3035(float(lon), float(lat)) ###NEEDS TO BE CHANGED
                    #print(crs,lon,lat)
                else:
                    crs_indicator= crs.replace("EPSG/0", "EPSG:")
                    lon,lat= trans4mfromwgs84(crs_indicator,float(lon), float(lat))
                #print(lon,lat)
                coverage_id = layer
                x1= float(lat)
                x2= x1 + resolution
                #print(x1,x2, resolution)
                #strlat="&subset=Lat({},{})"
                subset_lat = strlat.format(x1,x2)
                y1= float(lat)
                y2= y1 + resolution
                #strlon = "&subset=Long({},{})"
                subset_long= strlon.format(y1,y2)
                #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
                cov_id = "&COVERAGEID=" + coverage_id
                encode_format = "&FORMAT=text/csv"
                #response = requests.get(
                #base_wcs_url + "&request=GetCoverage" + cov_id + subset_long + subset_lat + encode_format,auth=(rasdaman_username, rasdaman_password), verify=False)
                #sample_result.append(response.content)
                print(layer,lat,lon, row["minlat"], row["maxlat"])
                if float(lat) > float(row["minlat"]) and float(lat) < float(row["maxlat"]) and float(lon) > float(row["minlong"]) and float(lon) < float(row["maxlong"] ):
                    my_dict[name_date]=set()
                    my_dict[name_date].add(layer)  ##this saves information on which
                else:
                    print("")


                ##add metadata link for description of the layer data

                ## remove some prints

                ## add to store values infile /make the response statement valid 

                

        #print("####")
