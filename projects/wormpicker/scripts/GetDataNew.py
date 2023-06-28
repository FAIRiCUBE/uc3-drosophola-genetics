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

my_dict={}
end_result=[]
fully_covered=[]
covered_spatially=[]
dropout=[]
with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/coveredsamples_wcs.csv", newline='') as csvfile:
    sample_reader = csv.DictReader(csvfile)
    for s_row in sample_reader:
        latitude_sample = s_row["lat"]
        longitude_sample = s_row["long"]
        name_date = s_row["sampleId"]
        date = str(name_date[-10:])
        time_asked= date + 'T00:00:00Z'
        #print("ORIGINAL COORDINATES OF SAMPLE\n:", name_date, latitude_original, longitude_original, "DATE:", time_asked)
        #datetime3 = parser.isoparse(time_asked)
        sample_result=[]
        #print(datetime3)
        layer=None
        converted_date=None
        my_dict[name_date]=[]  
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            search_time=parser.isoparse(converted_date)
            #print(search_time)
        except ValueError:
            converted_date = date
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
                latitude_request = latitude_sample
                longitude_request = longitude_sample
                time_min= row["f"]
                time_max= row["t"]
                #print("LIMITS: ", time_min, time_max, "MI\n", row["minlat"], row["maxlat"] )
                datetime1 = parser.isoparse(time_min)
                datetime2 = parser.isoparse(time_max)
                #print(layer, crs, lat, lon, resolution)
                #print(layer, search_time, converted_date, datetime1, datetime2) 
                #print(lat,lon)
                strlat="&subset=Y({},{})"
                strlon = "&subset=X({},{})"
                #print("SUBSETS DEFAULT",strlon,strlat)
                #print(abs(float(row["resolution_dim1"])))
                ansi_str="&subset=ansi(\"{}\")"
                if search_time > datetime1 and search_time < datetime2:
                    print("YES")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                    ##produce timestamp here
                    ansi_val=converted_date
                    subset_ansi=ansi_str.format(ansi_val)
                else:
                    print("")  # NOT COVERED TEMPORALLY
                    #set timestamp here a general one, the one were coverage is there
                    ansi_val=time_min
                    subset_ansi=ansi_str.format(ansi_val)     
                if abs(float(row["resolution_dim1"])) == abs(float(row["resolution_dim2"])):
                    resolution=abs(float(row["resolution_dim1"]))
                    #print("AT least axis resolution is identical")
                else:
                    #print("NOT same resolition in axis")
                    resolution=abs(float(row["resolution_dim2"])) ###still not sure why and what to do but spaceholder
                #print(resolution)
                #print(crs)
                if crs=="EPSG/0/4326":
                    latitude_request=latitude_request
                    longitude_request=longitude_request
                    strlat="&subset=Lat({},{})"
                    strlon = "&subset=Long({},{})"
                    print(crs,latitude_request,longitude_request)
                else:
                    crs_indicator= crs.replace("EPSG/0/", "EPSG:")
                    longitude_request,latitude_request= trans4mfromwgs84(crs_indicator,float(longitude_request), float(latitude_request))
                #print(strlat, strlon)
                #coverage_id = layer
                x1= float(latitude_request)
                x2= x1 + resolution
                #print(x1,x2, resolution)
                #strlat="&subset=Lat({},{})"
                subset_lat = strlat.format(x1,x2)
                y1= float(longitude_request)
                y2= y1 + resolution
                #strlon = "&subset=Long({},{})"
                subset_long= strlon.format(y1,y2)
                #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
                request_cov_id = "&COVERAGEID=" + layer
                request_encode_format = "&FORMAT=text/csv"
                #print(layer,lat,lon, row["minlat"], row["maxlat"])
                if float(latitude_request) > float(row["minlat"]) and float(latitude_request) < float(row["maxlat"]) and float(longitude_request) > float(row["minlong"]) and float(longitude_request) < float(row["maxlong"] ):
                    my_dict[name_date]=set()
                    my_dict[name_date].add(layer)  ##this saves information on which
                    response = requests.get(
                    base_wcs_url + "&request=GetCoverage" + request_cov_id + subset_ansi + subset_lat + subset_long + request_encode_format,auth=(rasdaman_username, rasdaman_password), verify=False)
                    sample_result.append((layer,response.text))
                else:
                    print("YAY")
                    dropout.append(("ERROR: Sample " + str(name_date) + " not covered in area of " + layer + ". You calculated EPSG for WGS84:" + str(latitude_request) + str(longitude_request) + ". Boundaries lat are" + row["minlat"]+ row["maxlat"] + "long: " + row["minlong"] + row["maxlong"]))
        print(sample_result)
        end_result.append((s_row["sampleId"], sample_result))
    print(end_result)


                ##add metadata link for description of the layer data

                ## remove some prints

                ## add to store values infile /make the response statement valid 

                

        #print("####")
