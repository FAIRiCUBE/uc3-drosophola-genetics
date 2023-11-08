import sys
import os 
import csv
import requests
from dotenv import dotenv_values
from dateutil import parser     
from datetime import datetime
#from module_crs_converter import trans4mEPSG
import math
from osgeo import gdal, osr


def trans4mEPSG(InputCRS,OutputCRS,y,x):
    src_crs = InputCRS
    src_srs = osr.SpatialReference()
    src_srs.SetFromUserInput(src_crs)
    tgt_crs = OutputCRS
    tgt_srs = osr.SpatialReference()
    tgt_srs.SetFromUserInput(tgt_crs)
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    x_t, y_t, z_t = transform.TransformPoint(x, y, 0)
    return x_t, y_t

#env_vars = dotenv_values()
## Access specific environment variables
#rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
#rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
#rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
#base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"

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
        filtered_data={}
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
                    filtered_data[row["sampleId"]] = sampleinfo
                    #filtered_data.append(sampleinfo)
        self.samples=filtered_data

##into README
#layer_info=[]
#cos=Coverage(layer_info)
#cos.getBoundary()
#cos.getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")

#samplescovered=cos.samples

#import os
#make this suitable for sample list
#def getSamples(path):
#    if os.path.exists(path):
#        if os.path.isfile(path):
#            with open(path, 'r') as file:
#                lines = file.readlines()
#                for line in lines:
#                    print(line.strip())  # Process each line as desired
#        else:
#            print(f"Error: '{path}' is not a file.")
#    else:
#        print(f"Error: Path '{path}' does not exist.")
#

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

#def show_selected_options():
#    selected_options = []
#    for index, var in enumerate(option_vars):
#        if var.get() == 1:
#            selected_options.append(index)
#    if selected_options:
#        for index in selected_options:
#            x = option_names[index]
#            messagebox.showinfo("Selected Option", f"You chose {x}")
#            chosen_layers.append(x)
#        window.destroy()
#mode="auto"
mode="manual"

#if mode=="manual":
#    objects_list = cos.layers
#    window = tk.Tk()
#    window.title("Select Objects")
#    option_vars = []
#    option_names=[]
#    chosen_layers=[]
#    for index, option in enumerate(objects_list):
#        var = tk.IntVar(window)
#        check_button = ttk.Checkbutton(window, text=option, variable=var)
#        check_button.grid(row=index, column=0, sticky="w")
#        option_vars.append(var)
#        option_names.append(option)
#    btn_submit = ttk.Button(window, text="Submit", command=show_selected_options)
#    btn_submit.grid(row=len(objects_list), column=0, pady=10)
#    window.mainloop()
#    use_layers=chosen_layers
#else:
#    use_layers=cos.layers
### now layers are in use_layers
#def show_selected_options():
#    selected_options = []
#    for index, var in enumerate(option_vars):
#        if var.get() == 1:
#            selected_options.append(index)
#    if selected_options:
#        for index in selected_options:
#            x = option_names[index]
#            messagebox.showinfo("Selected Option", f"You chose {x}")
#            chosen_layers.append(x)
#        window.destroy()

def select_objects(mode, object):
    cos=object
    if mode == "manual":
        def show_selected_options():
            global chosen_layers
            chosen_layers = [option_names[i] for i, var in enumerate(option_vars) if var.get() == 1]
            window.destroy()
        objects_list = cos.layers
        window = tk.Tk()
        window.title("Select Layers")
        option_vars = []
        option_names = []
        for index, option in enumerate(objects_list):
            var = tk.IntVar(window)
            check_button = ttk.Checkbutton(window, text=option, variable=var)
            check_button.grid(row=index, column=0, sticky="w")
            option_vars.append(var)
            option_names.append(option)
        btn_submit = ttk.Button(window, text="Submit", command=show_selected_options)
        btn_submit.grid(row=len(objects_list), column=0, pady=10)
        window.mainloop()
        use_layers = chosen_layers
    if mode == "automatic":
        use_layers = cos.layers
    return use_layers

###README USAGE
#layers_to_analyze=select_objects("automatic", cos)


##repeat this for samples!?!?!
###
#data_use=[]
#for layer in use_layers:
#    for i in range(0,len(cos._data)):
#        entries_list=list(cos._data[i])
#        try:
#            index=entries_list.index(layer)
#            data_use.append(entries_list)
#        except:
#            continue

def requestData(layerlist,samples, filepath):
    env_vars = dotenv_values()
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"
    result=[]
    header=[]
    header.append('Sample')
    for data_entry in range(0,len(layerlist)):
        layer=layerlist[data_entry][0]
        header.append(layer)
    log = open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/WormpickerRequest.log", "a")
    #sys.stdout = log
    for key, value in samples.items():
        sample_result=[]
        latitude_sample = value[0] 
        longitude_sample = value[1] 
        name_date = key
        date = str(name_date[-10:])
        time_asked= date + 'T00:00:00Z'
        sample_result.append(name_date)
        print(sample_result)
        #layer=None
        #converted_date=None
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            search_time=parser.isoparse(converted_date)
            print(search_time)
        except ValueError:
            converted_date = date
            continue
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry][0]
           #header.append(layer) 
            crs=layerlist[data_entry][1] 
            resolution=layerlist[data_entry][8] 
            latitude_request = latitude_sample
            longitude_request = longitude_sample
            time_min= layerlist[data_entry][6] 
            time_max= layerlist[data_entry][7]
            datetime1 = parser.isoparse(time_min)
            datetime2 = parser.isoparse(time_max)
            strlat="&subset=Y({},{})"
            strlon = "&subset=X({},{})"
            ansi_str="&subset=ansi(\"{}\")"    
            #print(latitude_request)
            if search_time > datetime1 and search_time < datetime2:
                print("YES")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                ##produce timestamp here
                ansi_val=converted_date
                subset_ansi=ansi_str.format(ansi_val)
            else:
                print(time_asked, "# NOT COVERED TEMPORALLY, only retrieving data for:", time_min)
                #set timestamp here a general one, the one were coverage is there
                ansi_val=time_min
                subset_ansi=ansi_str.format(ansi_val)
            if crs=="EPSG/0/4326":
                latitude_request=latitude_request
                longitude_request=longitude_request
                strlat="&subset=Lat({},{})"
                strlon = "&subset=Long({},{})"
            else:
                crs_indicator= crs.replace("EPSG/0/", "EPSG:")
                longitude_request,latitude_request= trans4mEPSG("EPSG:4326",crs_indicator,float(longitude_request), float(latitude_request))
                x1= float(latitude_request)
                x2= x1 + float(resolution)
                #x2=x1
                #print(x1,x2, resolution)
                #strlat="&subset=Lat({},{})"
                subset_lat = strlat.format(x1,x2)
                y1= float(longitude_request)
                #y2=y1
                y2= y1 + float(resolution)
                #strlon = "&subset=Long({},{})"
                subset_long= strlon.format(y1,y2)
                #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
            request_cov_id = "&COVERAGEID=" + layer
            request_encode_format = "&FORMAT=text/csv"
            #print(name_date, request_cov_id, request_encode_format, layer)
            if float(latitude_request) < float(layerlist[data_entry][3]) and float(latitude_request) > float(layerlist[data_entry][2]) and float(longitude_request) > float(layerlist[data_entry][4]) and float(longitude_request) < float(layerlist[data_entry][5]):
                response = requests.get(base_wcs_url + "&request=GetCoverage" + request_cov_id + subset_ansi + subset_lat + subset_long + request_encode_format,auth=(rasdaman_username, rasdaman_password), verify=False)
                valls=(str(response.text)[1:-1])
                sample_result.append(valls)
            else:
                #print("YAY")
                valls=("NA")
                sample_result.append(valls)
        result.append(sample_result)
    #sys.stdout = sys.__stdout__
    print("DONE")
    if filepath!="NONE":
        with open(filepath, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(result)
    else:
        return result

#
# x= requestData(data_use,cos.samples,"/adada")


#def getBoundary(data):
#    minlats=[]
#    maxlats=[]
#    minlongs=[]
#    maxlongs=[]  
#    for minl in range(1,len(data)):
#        crs=data[minl][1]
#        #print(minl)
#        if crs=="EPSG/0/4326":
#            min_y, min_x = trans4mEPSG("EPSG:4326","EPSG:3035", float(data[minl][4]), float(data[minl][2]))
#            max_y, max_x = trans4mEPSG("EPSG:4326","EPSG:3035",float(data[minl][5]), float(data[minl][3]))
#            minlats.append(min_x)
#            maxlats.append(max_x)
#            minlongs.append(min_y)
#            maxlongs.append(max_y)
#        else:
#            minlats.append(float(data[minl][2]))
#            #print(minlats)
#            maxlats.append(float(data[minl][3]))
#            #print(maxlats)
#            minlongs.append(float(data[minl][4]))
#            maxlongs.append(float(data[minl][5]))
#        minlat=min([x for x in minlats if x != float('inf')])
#        maxlat= max([x for x in maxlats if x != float('inf')])
#        minlong=min([x for x in minlongs if x != float('inf')])
#        maxlong=max([x for x in maxlongs if x != float('inf')])
#    return(minlat,maxlat,minlong,maxlong)
#
#s=getBoundary(layer_info)



###TO DO HERE:
# make it usable with new code

# geo extent in EPSG:4326
#xmin = 0
#xmax = 180
#ymin = 0
#ymax = 90
#
## geo resolution of Long (X) axis
#xres = 1
#
## geo resolution of Lat (Y) axis
## (e.g. Lat axis, then pixel 0th is at 90 degree and pixel 89th at 0 degree with yres = -1)
#yres = -1

def get_grid_indices_for_axis_X(geo_lower_bound, geo_upper_bound, xmin, xmax, xres):       
    print("geo_lower_bound_X: {}, geo_upper_bound_X: {}".format(geo_lower_bound, geo_upper_bound))
    lower_grid_index = math.floor( (geo_lower_bound - xmin) / xres )
    upper_grid_index = math.floor( (geo_upper_bound - xmin) / xres )
    if geo_upper_bound == xmax:
        upper_grid_index = upper_grid_index - 1       
    return "({}:{})".format(lower_grid_index, upper_grid_index)
    
def get_grid_indices_for_axis_Y(geo_lower_bound, geo_upper_bound, ymin, ymax, yres):
    print("geo_lower_bound_Y: {}, geo_upper_bound_Y: {}".format(geo_lower_bound, geo_upper_bound))
    lower_grid_index = math.floor( (geo_upper_bound - ymax) / yres )
    upper_grid_index = math.floor( (geo_lower_bound - ymax) / yres )
    if geo_lower_bound == ymin:
        upper_grid_index = upper_grid_index - 1   
    return "({}:{})".format(lower_grid_index, upper_grid_index)

def requestDataProcess(layerlist,samples, filepath):
    env_vars = dotenv_values()
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?service=WCS&version=2.1.0"
    result=[]
    header=[]
    header.append('Sample')
    for data_entry in range(0,len(layerlist)):
        layer=layerlist[data_entry][0]
        header.append(layer)
    log = open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/WormpickerRequest.log", "a")
    #sys.stdout = log
    for key, value in samples.items():
        sample_result=[]
        latitude_sample = value[0] 
        longitude_sample = value[1]
        name_date = key
        date = str(name_date[-10:])
        time_asked= date + 'T00:00:00Z'
        sample_result.append(name_date)
        print(sample_result)
        #layer=None
        #converted_date=None
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            search_time=parser.isoparse(converted_date)
            print(search_time)
        except ValueError:
            converted_date = date
            continue
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry][0]
           #header.append(layer) 
            crs=layerlist[data_entry][1] 
            latitude_request = latitude_sample
            longitude_request = longitude_sample
            time_min= layerlist[data_entry][6] 
            time_max= layerlist[data_entry][7]
            datetime1 = parser.isoparse(time_min)
            datetime2 = parser.isoparse(time_max)
            ymin = float(layerlist[data_entry][2]) 
            ymax = float(layerlist[data_entry][3]) 
            xmin = float(layerlist[data_entry][4]) 
            xmax = float(layerlist[data_entry][5]) 
            xres=float(layerlist[data_entry][8])
            yres=float(layerlist[data_entry][9])
            #strlat="&subset=Y({},{})"
            #strlon = "&subset=X({},{})"
            ansi_str="&subset=ansi(\"{}\")"    
            #print(latitude_request)
            if search_time > datetime1 and search_time < datetime2:
                print("YES")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                ##produce timestamp here
                ansi_val=converted_date
                #subset_ansi=ansi_str.format(ansi_val)
            else:
                print(time_asked, "# NOT COVERED TEMPORALLY, only retrieving data for:", time_min)
                #set timestamp here a general one, the one were coverage is there
                ansi_val=time_min
                #subset_ansi=ansi_str.format(ansi_val)
            if crs=="EPSG/0/4326":
                print("CRS IS 4326")
                y1=latitude_request
                x1=longitude_request
                #strlat="&subset=Lat({},{})"
                #strlon = "&subset=Long({},{})"
            else:
                print("CRS IS NOT WGS84")
                crs_indicator= crs.replace("EPSG/0/", "EPSG:")
                longitude_request,latitude_request= trans4mEPSG("EPSG:4326",crs_indicator,float(longitude_request), float(latitude_request))
                y1= float(latitude_request)
                #x2= x1 + float(resolution)
                #x2=x1
                #print(x1,x2, resolution)
                #strlat="&subset=Lat({},{})"
                #subset_lat = strlat.format(x1,x2)
                x1= float(longitude_request)
                #y2=y1
                #y2= y1 + float(resolution)
                #strlon = "&subset=Long({},{})"
                #subset_long= strlon.format(y1,y2)
                #print("CRS:",crs, "RESOLUTION:", resolution, subset_long, subset_lat)
                print(x1)
                grid_indices_axis_X = get_grid_indices_for_axis_X(x1, x1, xmin, xmax, xres)
                grid_indices_axis_Y = get_grid_indices_for_axis_Y(y1, y1, ymin, ymax, yres)
                output_format = "text/csv"
            #print(name_date, request_cov_id, request_encode_format, layer)
            if x1 < xmax and x1 > xmin and y1 > ymin and y1 < ymax:
                print("YES CONDITIONS ARE MET")
                query = f"for $c in ({layer}) return encode($c[ ansi(\"{ansi_val}\"), X:\"CRS:1\"{grid_indices_axis_X}, Y:\"CRS:1\"{grid_indices_axis_Y} ], \"{output_format}\")"
                # Construct the URL with variables
                url = f"{base_wcs_url}&REQUEST=ProcessCoverages&QUERY={query}"
                response = requests.get(url, auth=(rasdaman_username, rasdaman_password))
                valls=(str(response.text)[1:-1])
                sample_result.append(valls)
                print(query)
            else:
                print("NAY")
                valls=("NA")
                sample_result.append(valls)
        result.append(sample_result)
    #sys.stdout = sys.__stdout__
    print("----------------------")
    print("                      ")
    if filepath!="NONE":
        with open(filepath, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(result)
    else:
        return result
