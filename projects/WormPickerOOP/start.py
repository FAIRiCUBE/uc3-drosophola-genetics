
from module_crs_converter import trans4mEPSG
import csv 

class Coverage(object):
    def __init__(self,data):
        self._data=data
        self._ID= data[0]
        self._CRS=data[1]
        self._minlat=data[2]
        self._maxlar= data[3]
        self._minlong=data[4]
        self._maxlong=data[5]
        self._fromtime=data[6] 
        self._totime=data[7]
        lnames=[]
        for lna in range(1,len(data)):
            lnames.append(data[lna][0])
        self.layers=lnames
    def getBoundary(self):
        minlats=[]
        maxlats=[]
        minlongs=[]
        maxlongs=[]  
        for minl in range(1,len(self._data)):
            crs=self._data[minl][1]
            if crs=="EPSG/0/4326":
                min_y, min_x = trans4mEPSG("EPSG:4326","EPSG:3035", float(self._data[minl][4]), float(self._data[minl][2]))
                max_y, max_x = trans4mEPSG("EPSG:4326","EPSG:3035",float(self._data[minl][5]), float(self._data[minl][3]))
                minlats.append(min_x)
                maxlats.append(max_x)
                minlongs.append(min_y)
                maxlongs.append(max_y)
            else:
                minlats.append(float(self._data[minl][2]))
                maxlats.append(float(self._data[minl][3]))
                minlongs.append(float(self._data[minl][4]))
                maxlongs.append(float(self._data[minl][5]))
        self.minlat=min([x for x in minlats if x != float('inf')])
        self.maxlat= max([x for x in maxlats if x != float('inf')])
        self.minlong=min([x for x in minlongs if x != float('inf')])
        self.maxlong=max([x for x in maxlongs if x != float('inf')])
    def getSamples(self,path):
        filtered_data=[]
        with open(path, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                lat=row["lat"]
                long=row["long"]
                if lat == 'NA' and long == 'NA':
                    continue
                else:
                    long,lat=trans4mEPSG("EPSG:4326","EPSG:3035",float(long),float(lat))
                if float(lat) > self.minlat and float(lat) < self.maxlat and float(long) > self.minlong and float(long) < self.maxlong:
                    sampleinfo=(row["lat"],row["long"])
                    filtered_data.append(row["sampleId"])
                    filtered_data.append(sampleinfo)
        self.samples=filtered_data


cos=Coverage(layer_info)
cos.getBoundary()
cos.getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")

#samplescovered=cos.samples

import os
#make this suitable for sample list
def getSamples(path):
    if os.path.exists(path):
        if os.path.isfile(path):
            with open(path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    print(line.strip())  # Process each line as desired
        else:
            print(f"Error: '{path}' is not a file.")
    else:
        print(f"Error: Path '{path}' does not exist.")


#getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")

#def getCoverageData():
    ###request and store??


## 1. Get Layers (shall download the layerinfo from the requested server)
## 2. Get Samples (give an input file to create list/object of samples)
## 3. Filter
## 4. GetCoverage (is a function that combines download)
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox



def show_selected_options():
    selected_options = []
    for index, var in enumerate(option_vars):
        if var.get() == 1:
            selected_options.append(index)
    if selected_options:
        for index in selected_options:
            x = option_names[index]
            messagebox.showinfo("Selected Option", f"You chose {x}")
            chosen_layers.append(x)
        window.destroy()


#cos=Coverage(layer_info)
#objects_list = ['AverageChlorophyll', 'corine_land_cover', 'dominant_leaf_type', 'european_settlement_map', 'grassland__change', 'grassland_status', 'grassland_status_demo', 'high_resolution_layer_forest_type', 'imperviousness', 'maxes_sentinel2_2018_flevopolder_10m_7x4bands', 'near_surface_meteorological_variables_ASurf', 'sentinel2_2018_flevopolder_10m_7x4bands', 'small_woody_features_16', 'small_woody_features_32', 'small_woody_features_4', 'small_woody_features_64', 'small_woody_features_8', 'tree_cover_density', 'water_and_wetness']
objects_list = cos.layers
window = tk.Tk()
window.title("Select Objects")

option_vars = []
option_names=[]
chosen_layers=[]
for index, option in enumerate(objects_list):
    var = tk.IntVar()
    check_button = ttk.Checkbutton(window, text=option, variable=var)
    check_button.grid(row=index, column=0, sticky="w")
    option_vars.append(var)
    option_names.append(option)

btn_submit = ttk.Button(window, text="Submit", command=show_selected_options)
btn_submit.grid(row=len(objects_list), column=0, pady=10)

window.mainloop()

### now layers are in chosen layers
chosen_layers


#def requestData(layerlist,samplelist):
#    for

for layer in chosen_layers:
    for i in range(0,len(cos._data)):
        entries_list=list(cos._data[i])
        try:
            index=entries_list.index(layer)
            #print(i,index)
            ###insert###
            f any(value == "NA" for value in row.values()):
                    continue
                layer=row["CoverageID"]
                #layerlist.append(layer)
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
                #print(layer, crs,resolution)
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
                    #print(crs,longitude_request,latitude_request)
                else:
                    crs_indicator= crs.replace("EPSG/0/", "EPSG:")
                    longitude_request,latitude_request= trans4mEPSG("EPSG:4326",crs_indicator,float(longitude_request), float(latitude_request))
                    #print(longitude_request,latitude_request)
                #print(strlat, strlon)
                #coverage_id = layer
                x1= float(latitude_request)
                #x2= x1 + resolution
                x2=x1
                #print(x1,x2, resolution)
                #strlat="&subset=Lat({},{})"
                subset_lat = strlat.format(x1,x2)
                y1= float(longitude_request)
                y2=y1
                #y2= y1 + resolution
                #strlon = "&subset=Long({},{})"
                subset_long= strlon.format(y1,y2)
                #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
                request_cov_id = "&COVERAGEID=" + layer
                request_encode_format = "&FORMAT=text/csv"
                if float(latitude_request) > float(row["minlat"]) and float(latitude_request) < float(row["maxlat"]) and float(longitude_request) > float(row["minlong"]) and float(longitude_request) < float(row["maxlong"] ):
                    my_dict[name_date]=set()
                    my_dict[name_date].add(layer)  ##this saves information on which
                    response = requests.get(
                    base_wcs_url + "&request=GetCoverage" + request_cov_id + subset_ansi + subset_lat + subset_long + request_encode_format,auth=(rasdaman_username, rasdaman_password), verify=False)
                    valls=(str(response.text)[1:-1])
                    sample_result.append(valls)
                else:
                    #print("YAY")
                    valls=("NA")
                    sample_result.append(valls)
                    
            ####
        except ValueError:
            continue