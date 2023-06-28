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
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

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


### 2. set a starting point  

### Testing 0/0 corner of grid cell 

T1=(900525,1400525)
T2=(900525,1400550) 
T3=(900525,1400575) 
T4=(900550,1400525)
T5=(900550,1400550) 
T6=(900550,1400575) 
T7=(900575,1400525) 
T8=(900575,1400550) 
T9=(900575,1400575)

 ### Q2 at 0.5
layer=getlayer()

points=["T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9"] 

#points=["T5", "T6", "T7", "T8", "T9"] 

for Testpoint in points:
    #print(Testpoint)
    coords=globals()[Testpoint]
    #print(coords)
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
    ### 3. define directions that will be testet
    Qs=((1,1),(-1,1),(-1,-1),(1,-1))
    ### 4. test all and get results
    facts=[]
    values=[]
    resolutions=[]
    cellnumbers=[]
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
    ### 3. define directions that will be testet
    Qs=((1,1),(-1,1),(-1,-1),(1,-1))
    ### 4. test all and get results
    facts=[]
    values=[]
    resolutions=[]
    cellnumbers=[]
    cellnumbersquadr=[]
    exes=[] 
    for k in range (0,100):
        for i in range(0,4):
            x1= float(longitude_request)
            x2= x1 + resolution*Qs[i][0]*k*0.25
            y1= float(latitude_request)
            y2= y1 + resolution*Qs[i][1]*k*0.25
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
            facts.append(("Quadrant:",i,"Resolution", k*0.25, "Exclam:", len(xi), "Within:", cells, "GridBoundary:", y1,y2, x1,x2))
            resolutions.append(k*0.25)
            cellnumbers.append(cells)
            exes.append((len(xi),cells))
            values.append((k*0.25, i, value))
        exes.append(k*0.25)
            #print(value)
    ### 5. table results 
    ### 6. find pattern in result 
    ########## REALLLLY
    for i in range(0, len(facts), 4):
    #for i in range(0, 20, 4):
        q1 = i
        q2 = i + 1
        q3 = i + 2
        q4 = i + 3
        resolutionf = facts[i][3]
        x1, y1 = facts[q1][5], facts[q1][7]
        x2, y2 = facts[q2][5], facts[q2][7]
        x3, y3 = facts[q3][5], facts[q3][7]
        x4, y4 = facts[q4][5], facts[q4][7]
        bb1=facts[i][9]
        bb2=facts[i][10]
        bb3=facts[i][11]
        bb4=facts[i][12]
        # Create a new figure and axis
        fig, ax = plt.subplots()
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        # Create FancyArrowPatch objects for the arrows
        arrow1 = FancyArrowPatch((0, 0), (resolutionf, 0), color='red', arrowstyle='->', linewidth=2)
        arrow2 = FancyArrowPatch((0, 0), (0, resolutionf), color='red', arrowstyle='->', linewidth=2)
        arrow3 = FancyArrowPatch((0, 0), (-resolutionf, 0), color='red', arrowstyle='->', linewidth=2)
        arrow4 = FancyArrowPatch((0, 0), (0, -resolutionf), color='red', arrowstyle='->', linewidth=2)
        # Add the arrows to the plot
        ax.add_patch(arrow1)
        ax.add_patch(arrow2)
        ax.add_patch(arrow3)
        ax.add_patch(arrow4)
        # Draw the squares as rectangles
        square1 = Rectangle((0, 0), x1, y1, edgecolor='blue',  facecolor=(0, 0, 1, 0.5))
        square2 = Rectangle((0, 0), -x2, y2, edgecolor='red', facecolor='none')
        square3 = Rectangle((0, 0), -x3, -y3, edgecolor='green', facecolor='none')
        square4 = Rectangle((0, 0), x4, -y4, edgecolor='magenta', facecolor='none')
        # Add the rectangles to the plot
        ax.add_patch(square1)
        ax.add_patch(square2)
        ax.add_patch(square3)
        ax.add_patch(square4)
        ax.set_aspect('equal')
        ax.grid(visible=True, linestyle=':')
        ax.set_xlim(-(resolutionf+5), resolutionf+5)
        ax.set_ylim(-(resolutionf+5), resolutionf+5)
        # Add a legend
        legend_labels =[f'Resolution factor={resolutionf} for {resolution} metre ', 'Q1', 'Q2', 'Q3', 'Q4', f'Coordinate Bounding ={bb1}-{bb2}:{bb3}-{bb4}']
        ax.legend([arrow1, square1, square2, square3, square4, square1], legend_labels)
        # Display all the plots
        save_path = f'/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/EggTest/output/_new_{Testpoint}_{resolutionf}.png'
        plt.savefig(save_path)    
        # Show the plot
        #plt.show()


####
cellnumbersquadr=[]
exes=[] 
for k in range (0,100):
    for i in range(0,4):
        x1= float(longitude_request)
        x2= x1 + resolution*Qs[i][0]*k*0.25
        y1= float(latitude_request)
        y2= y1 + resolution*Qs[i][1]*k*0.25
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
        facts.append(("Quadrant:",i,"Resolution", k*0.25, "Exclam:", len(xi), "Within:", cells, "GridBoundary:", y1,y2, x1,x2))
        resolutions.append(k*25)
        cellnumbers.append(cells)
        exes.append((len(xi),cells))
        values.append((k*0.25, i, value))
    exes.append(k*0.25)
        #print(value)




### 5. table results 

### 6. find pattern in result 

########## REALLLL

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

for i in range(0, len(facts), 4):
#for i in range(0, 40, 4):
    q1 = i
    q2 = i + 1
    q3 = i + 2
    q4 = i + 3
    resolutionf = facts[i][3]
    x1, y1 = facts[q1][5], facts[q1][7]
    x2, y2 = facts[q2][5], facts[q2][7]
    x3, y3 = facts[q3][5], facts[q3][7]
    x4, y4 = facts[q4][5], facts[q4][7]
    bb1=facts[i][9]
    bb2=facts[i][10]
    bb3=facts[i][11]
    bb4=facts[i][12]
    # Create a new figure and axis
    fig, ax = plt.subplots()
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    # Create FancyArrowPatch objects for the arrows
    arrow1 = FancyArrowPatch((0, 0), (resolutionf, 0), color='red', arrowstyle='->', linewidth=2)
    arrow2 = FancyArrowPatch((0, 0), (0, resolutionf), color='red', arrowstyle='->', linewidth=2)
    arrow3 = FancyArrowPatch((0, 0), (-resolutionf, 0), color='red', arrowstyle='->', linewidth=2)
    arrow4 = FancyArrowPatch((0, 0), (0, -resolutionf), color='red', arrowstyle='->', linewidth=2)
    # Add the arrows to the plot
    ax.add_patch(arrow1)
    ax.add_patch(arrow2)
    ax.add_patch(arrow3)
    ax.add_patch(arrow4)
    # Draw the squares as rectangles
    square1 = Rectangle((0, 0), x1, y1, edgecolor='blue',  facecolor=(0, 0, 1, 0.5))
    square2 = Rectangle((0, 0), -x2, y2, edgecolor='red', facecolor='none')
    square3 = Rectangle((0, 0), -x3, -y3, edgecolor='green', facecolor='none')
    square4 = Rectangle((0, 0), x4, -y4, edgecolor='magenta', facecolor='none')
    # Add the rectangles to the plot
    ax.add_patch(square1)
    ax.add_patch(square2)
    ax.add_patch(square3)
    ax.add_patch(square4)
    ax.set_aspect('equal')
    ax.grid(visible=True, linestyle=':')
    ax.set_xlim(-(resolutionf+5), resolutionf+5)
    ax.set_ylim(-(resolutionf+5), resolutionf+5)
    # Add a legend
    legend_labels =[f'Resolution factor={resolutionf} for {resolution} metre ', 'Q1', 'Q2', 'Q3', 'Q4', f'Coordinate Bounding ={bb1}-{bb2}:{bb3}-{bb4}']
    ax.legend([arrow1, square1, square2, square3, square4, square1], legend_labels)
    # Display all the plots
    save_path = f'/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/EggTest/output/Fake_k1_{resolutionf}.png'
    plt.savefig(save_path)    
    # Show the plot
    #plt.show()



#redo for other alyer and other sample point 

r_set1=resolutions
c_set1=cellnumbers
plt.scatter(r_set1, c_set1, color='blue', marker='o', label='Austria')
r_set2=resolutions
c_set2=cellnumbers
plt.scatter(r_set2, c_set2, color='red', marker='o', label='Spain')
r_set3=resolutions
c_set3=cellnumbers
plt.scatter(r_set3, c_set3, color='red', marker='o', label='Serbia')
plt.show()






