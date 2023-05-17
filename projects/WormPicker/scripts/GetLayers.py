import requests
from dotenv import dotenv_values
import json
import xmltodict

layer_info = []
layer_info.append(("CoverageID", "CRS", "minlat","maxlat", "minlong", "maxlong","f","t"))
# Load environment variables from the .env file
env_vars = dotenv_values()

# Access specific environment variables
rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')

base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"
response = requests.get(base_wcs_url + "&request=GetCapabilities", auth=(rasdaman_username, rasdaman_password))

wcs_capabilities = xmltodict.parse(response.content)
wcs_coverage_summary = wcs_capabilities['wcs:Capabilities']['wcs:Contents']['wcs:CoverageSummary']
#print(json.dumps(wcs_coverage_summary, indent=2))
type(wcs_coverage_summary)

for i in range(0,len(wcs_coverage_summary)):
    coverage_id= wcs_coverage_summary[i]['wcs:CoverageId']
    dim=int(wcs_coverage_summary[i]['ows:BoundingBox']['@dimensions'])
    crs_str = wcs_coverage_summary[i]['ows:BoundingBox']['@crs'][-11:]
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
        layer_info.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_min))
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

with open("output/layer_info_WCS.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    #writer.writerow(["LayerName", "Min X", "Max X", "Min Y", "Max Y", "TimeStart", "TimeEnd"])
    writer.writerows(layer_info)
