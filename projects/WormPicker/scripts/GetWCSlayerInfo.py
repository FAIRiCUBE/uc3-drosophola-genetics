###try

import re
import requests
import xml.etree.ElementTree as ET
import csv 


#url = "https://ows.rasdaman.org/rasdaman/ows?service=WCS&request=GetCapabilities"

from dotenv import dotenv_values

# Load environment variables from the .env file
env_vars = dotenv_values()

# Access specific environment variables
rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')

base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"
response = requests.get(base_wcs_url + "&request=GetCapabilities", auth=(rasdaman_username, rasdaman_password))

# Send the GetCapabilities request and get the response
#response = requests.get(url)

# Parse the response as an XML document
tree = ET.fromstring(response.content)

# Find all CoverageId elements in the XML document

coverages = tree.findall(".//{http://www.opengis.net/wcs/2.1}CoverageSummary")

layer_info = []

## BRAND NEW  
# Loop through the CoverageId elements and print the ones that match the desired namespaces
for layer in coverages:
    coverage_id=layer.findall(".//{http://www.opengis.net/wcs/2.1}CoverageId")
    for title in coverage_id:
        #print(title.text)
        layer_name=title.text
    bbwg48 = layer.findall(".//{http://www.opengis.net/ows/2.0}WGS84BoundingBox")
    for x in bbwg48:
        min= x.findall(".//{http://www.opengis.net/ows/2.0}LowerCorner")
        for y in min:
        #print(y.tag)
            minlat= y.text.split(" ")[1]
            minlong=y.text.split(" ")[0]
        max= x.findall(".//{http://www.opengis.net/ows/2.0}UpperCorner")
        for z in max:
        #print(y.tag)
            maxlat= z.text.split(" ")[1]
            maxlong=z.text.split(" ")[0]
    #print(minlat, maxlat, minlong, maxlong)
    #layer_info.append((minlat,maxlat, minlong, maxlong))
    bb = layer.findall(".//{http://www.opengis.net/ows/2.0}BoundingBox")
    for x in bb:
        stamp=x.attrib['crs'].split('&')[0]
    #match=re.search(r'(?<=\?1=).*AnsiDate', stamp)
    if (x.attrib.get('dimensions') in ['3', '4']) and re.search(r'(?<=\?1=).*AnsiDate', x.attrib['crs']):
        lc = x.findall(".//{http://www.opengis.net/ows/2.0}LowerCorner")
        for i in lc: 
            #print("")
            #print(v.text)
            if re.search(r'(?<=\?1=).*AnsiDate', stamp):
                f=i.text.split(" ")[0] ## timestamp
            else:
                f=i.text.split(" ")[2] ## timestamp
        uc = x.findall(".//{http://www.opengis.net/ows/2.0}UpperCorner")
        for v in uc: 
            #print("")
            if re.search(r'(?<=\?1=).*AnsiDate', stamp):
                t=v.text.split(" ")[0] ## timestamp
            else:
                t=v.text.split(" ")[2] ## timestamp
        #print(f,t)
        #layer_info.append(f)
        #layer_info.append(t)
        print("1")
    else:
        t,f=("NA","NA")
        #layer_info.append('"NA" "NA"')
    layer_info.append((layer_name, minlat,maxlat, minlong, maxlong,f,t))


with open("output/fc_layer_info_WCS.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["LayerName", "Min X", "Max X", "Min Y", "Max Y", "TimeStart", "TimeEnd"])
    writer.writerows(layer_info)





