import requests
from dotenv import dotenv_values
import json
import xmltodict
import os
import csv


#create an empty array and add headers
layer_info = []
layer_info.append(("CoverageID", "CRS", "minlat","maxlat", "minlong", "maxlong","f","t", "resolution_dim1", "resolution_dim2", "resolution_dim3/time"))

#Load environment variables from the .env file
env_vars = dotenv_values()

# Access specific environment variables (useful when declared in a unix environment already)
rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')

base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"
response = requests.get(base_wcs_url + "&request=GetCapabilities", auth=(rasdaman_username, rasdaman_password))

wcs_capabilities = xmltodict.parse(response.content)
wcs_coverage_summary = wcs_capabilities['wcs:Capabilities']['wcs:Contents']['wcs:CoverageSummary']
#print(json.dumps(wcs_coverage_summary, indent=2))
type(wcs_coverage_summary)

cov_id_list=[]

#extract all coverage_ids and informations from the service endpoint (working on getCapabilities results)
for i in range(0,len(wcs_coverage_summary)):
    coverage_id= wcs_coverage_summary[i]['wcs:CoverageId']
    dim=int(wcs_coverage_summary[i]['ows:BoundingBox']['@dimensions'])
    #crs_str = wcs_coverage_summary[i]['ows:BoundingBox']['@crs'][-11:]
    crs_str = wcs_coverage_summary[i]['ows:BoundingBox']['@crs']
    crs = '/'.join(crs_str.split('/')[-3:])
    if dim==3:
        bb_low=wcs_coverage_summary[i]['ows:BoundingBox']['ows:LowerCorner']
        bb_upp=wcs_coverage_summary[i]['ows:BoundingBox']['ows:UpperCorner']
        x_min=bb_low.split(" ")[1] 
        y_min=bb_low.split(" ")[2]
        date_min=bb_low.strip().split(" ")[0][1:-1]
        date_max=bb_upp.strip().split(" ")[0][1:-1]
        x_max=bb_upp.split(" ")[1] 
        y_max=bb_upp.split(" ")[2]
        #print(coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_min)
        layer_info.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_max))
    else:
        bb_low=wcs_coverage_summary[i]['ows:BoundingBox']['ows:LowerCorner']
        bb_upp=wcs_coverage_summary[i]['ows:BoundingBox']['ows:UpperCorner']
        x_min=bb_low.split(" ")[0] 
        y_min=bb_low.split(" ")[-1+dim]
        x_max=bb_upp.split(" ")[0] 
        y_max=bb_upp.split(" ")[-1+dim]
        date_min="NA"
        date_max="NA"
        #print(coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_min)
        layer_info.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_min))

##add resolution via describe coverage

for ID in range(1,len(layer_info)):
#for ID in range(3,4):
    coverage=layer_info[ID][0]
    #print(coverage)
    response = requests.get(base_wcs_url + "&request=DescribeCoverage&coverageId=" + coverage,
                        auth=(rasdaman_username, rasdaman_password)
                        )
    wcs_coverage_description = xmltodict.parse(response.content)
    #print(json.dumps(wcs_coverage_description, indent=2))
    rr=[]
    try:
        l=len(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gml:RectifiedGrid']['gml:offsetVector'])
        for i in reversed(range(0,l)):
            resi=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gml:RectifiedGrid']['gml:offsetVector'][i]['#text']
            resolution=resi.split(" ")[i]
            rr.append(resolution)
            #x=tuple([resolution])
        layer_info[ID]=layer_info[ID]+tuple(rr)
        #print(layer_info[ID])
    except KeyError:
        pass
    try:
        l=len(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gmlrgrid:ReferenceableGridByVectors']['gmlrgrid:generalGridAxis'][0]['gmlrgrid:GeneralGridAxis']['gmlrgrid:offsetVector']['#text'].split(" "))
        for i in reversed(range(0,l)):
            resi=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gmlrgrid:ReferenceableGridByVectors']['gmlrgrid:generalGridAxis'][i]['gmlrgrid:GeneralGridAxis']['gmlrgrid:offsetVector']['#text'] 
            resolution=resi.split(" ")[i]
            rr.append(resolution)
            #x=tuple([resolution])
            #print(x)
        layer_info[ID]=layer_info[ID]+tuple(rr)
        #print(layer_info[ID])
    except KeyError:
        pass
    try:
        wcs_coverage_description['ows:ExceptionReport']
        rr=["NA", "NA", "NA"]
        layer_info[ID]=layer_info[ID]+tuple(rr)
    except KeyError:
        pass

os.chdir("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output")
with open("layer_info_WCS.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(layer_info)


