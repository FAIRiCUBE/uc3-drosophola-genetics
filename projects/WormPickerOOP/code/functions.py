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

env_path = '/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/WormPickerOOP'

# Load environment variables from the file
with open(env_path) as f:
    for line in f:
        key, value = line.strip().split('=', 1)
        os.environ[key] = value

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

def select_objects(mode, object, items_per_page=15):
    cos = object
    if mode == "manual":
        def show_selected_options():
            nonlocal chosen_layers
            chosen_layers = [option_names[i] for i, var in enumerate(option_vars) if var.get() == 1]
            window.destroy()
        def next_page():
            nonlocal current_page
            current_page += 1
            show_page()
        def prev_page():
            nonlocal current_page
            if current_page > 0:
                current_page -= 1
                show_page()
        def show_page():
            nonlocal next_button, prev_button, btn_submit
            for widget in window.winfo_children():
                widget.grid_forget()
            start = current_page * items_per_page
            end = (current_page + 1) * items_per_page
            for index, option in enumerate(objects_list[start:end]):
                var = tk.IntVar(value=0)
                check_button = ttk.Checkbutton(window, text=option, variable=var)
                check_button.grid(row=index, column=0, sticky="w")
                option_vars.append(var)
                option_names.append(option)
            if start > 0:
                prev_button.grid(row=items_per_page + 1, column=0, pady=10, padx=5, sticky="ew")
            if end < len(objects_list):
                next_button.grid(row=items_per_page + 1, column=1, pady=10, padx=5, sticky="ew")
            btn_submit.grid(row=items_per_page + 2, column=0, pady=10, sticky="ew")
        objects_list = cos.layers
        chosen_layers = []
        current_page = 0
        window = tk.Tk()
        window.title("Select Layers")
        window.geometry("800x400")  # Set the size of the window
        window.resizable(False, False)  # Make the window non-resizable
        option_vars = []
        option_names = []
        next_button = ttk.Button(window, text="Next", command=next_page)
        prev_button = ttk.Button(window, text="Previous", command=prev_page)
        btn_submit = ttk.Button(window, text="Submit", command=show_selected_options)
        show_page()
        window.mainloop()
        use_layers = chosen_layers
    elif mode == "automatic":
        use_layers = cos.layers
    else:
        raise ValueError("Invalid mode provided.")
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

def requestData(layerlist,samples,filepath):
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
    log = open("WormpickerRequest.log", "a")
    sys.stdout = log
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
    if filepath !="NONE":
        f=filepath
        with open(f, "w", newline="") as csvfile:
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

def is_valid_iso_date(date_string):
    try:
        # Attempt to parse the string as an ISO-formatted date
        datetime.fromisoformat(date_string)
        return True
    except ValueError:
        return False

def process_date_string(date_string):
    if is_valid_iso_date(date_string):
        print(f"{date_string} is already in ISO format. Doing nothing.")
        return date_string
    else:
        # Add your processing logic here if the date_string is not in ISO format
        print(f"Processing {date_string}...")
        # For now, just returning the input date_string as an example
        return date_string

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


def round_down_to_day(date_str):
    # Parse the date string
    try:
        date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        try:
            date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H:%M")
        except ValueError:
            date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H")
    rounded_down_date = date_obj.replace(hour=0, minute=0, second=0)
    rounded_down_date_str = rounded_down_date.strftime("%Y-%m-%dT%H:%M:%S")
    return rounded_down_date_str

def requestDataProcess(infoheader,layerlist,samples, filepath):
    env_vars = dotenv_values()
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?SERVICE=WCS&VERSION=2.0.1"
    result=[]
    header=[]
    header.append('Sample')
    if len(layerlist) > 1:
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry][0]
            header.append(layer)
    else:
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry]
            header.append(layer)
    log=open("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/WormpickerRequest.log", "a")
    sys.stdout = log
    all_results=[] 
    for key, value in samples.items():
        sample_result=[]
        latitude_sample=value[0] 
        longitude_sample=value[1]
        name_date = key
        print(" ")
        print(" ")
        print("NEW SAMPLE")
        print("COORDINATES:", latitude_sample, longitude_sample)
        try:
            date = value[2]
            print("DATE ASKED:",date)
            type(date)
        except IndexError:
            date = str(name_date[-10:])
            print("DATE ASKED:",date)
        try:
            daydate=round_down_to_day(date)
            converted_date=datetime.strptime(daydate, "%Y-%m-%dT%H:%M:%S").isoformat() + "Z"
            time_asked=converted_date
            #print("correct format")
        except ValueError:
            date = process_date_string(date + 'T00:00:00Z')
            time_asked=str(date)
            print(time_asked, name_date)
            #print(sample_result)
        sample_result.append(name_date)
        sample_result.append(time_asked)
        all_results.append(sample_result)
        #print(sample_result)
        print("QUERYING:", time_asked)
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            search_time=parser.isoparse(converted_date)
            print(search_time)
        except ValueError:
            print("converted_cate=date")
            search_time=parser.isoparse(time_asked)
        for data_entry in range(0,len(layerlist)):
            i_header=infoheader
            print(i_header)
            latitude_request = latitude_sample
            longitude_request = longitude_sample
            layer=layerlist[data_entry][i_header.index("CoverageID")]
            crs=layerlist[data_entry][i_header.index("CRS")]
            time_min= layerlist[data_entry][i_header.index("f")]
            time_max= layerlist[data_entry][i_header.index("t")]
            time_min=time_min.replace('"','')
            time_max=time_max.replace('"','')
            datetime1 = parser.isoparse(time_min)
            datetime2 = parser.isoparse(time_max)
            AxisLabels=layerlist[data_entry][i_header.index("axislabels")]
            TimeLabel=AxisLabels.split(",")[0]
            YLabel=AxisLabels.split(",")[1]
            XLabel=AxisLabels.split(",")[2]
            print(TimeLabel,YLabel, XLabel)
            ymin= float(layerlist[data_entry][i_header.index("minlat")])
            ymax = float(layerlist[data_entry][i_header.index("maxlat")])
            xmin= float(layerlist[data_entry][i_header.index("minlong")])
            xmax = float(layerlist[data_entry][i_header.index("maxlong")])
            xres=float(layerlist[data_entry][i_header.index("resolution_dim1")])
            yres=float(layerlist[data_entry][i_header.index("resolution_dim2")])
            #strlat="&subset=Y({},{})"
            #strlon = "&subset=X({},{})"
            ansi_str_o="&subset=ansi(\"{}\")"
            ansi_str_n=ansi_str_o.replace("ansi",TimeLabel)
            if search_time > datetime1 and search_time < datetime2:
                print(search_time, "- COVERED TEMPROALLY")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                ##produce timestamp here
                ansi_val=time_asked
            else:
                min_date = datetime1 if abs((datetime1 - search_time).days) < abs((datetime2 - search_time).days) else datetime2
                dist=abs(min_date-search_time).days
                print(time_asked, "NOT COVERED TEMPORALLY, querying:", min_date, ". This is", dist, "days apart!")
                #set timestamp here a general one, the one were coverage is there
                ansi_val=time_min
            if crs=="EPSG/0/4326":
                print("CRS:", crs)
                y1=float(latitude_request)
                x1=float(longitude_request)
                #strlat="&subset=Lat({},{})"
                #strlon = "&subset=Long({},{})"
            else:
                print(layer)
                crs_indicator= crs.replace("EPSG/0/", "EPSG:")
                #print("CRS IS NOT WGS84. It is:", crs_indicator)
                longitude_request,latitude_request= trans4mEPSG("EPSG:4326",crs_indicator,float(longitude_request), float(latitude_request))
                y1= float(latitude_request)
                x1= float(longitude_request)
            grid_indices_axis_X = get_grid_indices_for_axis_X(x1, x1, xmin, xmax, xres)
            grid_indices_axis_Y = get_grid_indices_for_axis_Y(y1, y1, ymin, ymax, yres)
            print("Calulculating Grid Indices from Coordinates:", x1, y1)
            print(grid_indices_axis_X, xres)
            print(grid_indices_axis_Y, yres)
            #corX= abs(int(grid_indices_axis_X[:-1][1:].split(":")[0]))
            #corY= abs(int(grid_indices_axis_Y[:-1][1:].split(":")[0]))
            #print(corX,corX)
            #print(corY,corY)
            grIx=("("+str(abs(int(grid_indices_axis_X[:-1][1:].split(":")[0]))) + ":" + str(abs(int(grid_indices_axis_X[:-1][1:].split(":")[0]))))+")"
            grIy=("("+str(abs(int(grid_indices_axis_Y[:-1][1:].split(":")[0]))) + ":" + str(abs(int(grid_indices_axis_Y[:-1][1:].split(":")[0]))))+")"
            output_format = "text/csv"
            print("Corrected Grid Indices:", grIx, grIy)
            #print(name_date, request_cov_id, request_encode_format, layer)
            if x1 <= xmax and x1 >= xmin and y1 >= ymin and y1 <= ymax and ansi_val != time_min and ansi_val != time_max:
                print("The Layer:", layer, " is within temporal coverage. Querying for", time_asked, ansi_val)
                #TimeLabel="date"
                #XLabel="Lon"
                #YLabel="Lat"
                query = f"for $c in ({layer}) return encode($c[{TimeLabel}(\"{ansi_val}\"),{XLabel} :\"CRS:1\"{grIx},{YLabel} :\"CRS:1\"{grIy} ], \"{output_format}\")"
                #query = f"for $c in ({layer}) return encode($c[date:(\"{ansi_val}\"),Lon :\"CRS:1\"{grid_indices_axis_X},Lat :\"CRS:1\"{grid_indices_axis_Y} ], \"{output_format}\")"
                # Construct the URL with variables
                url = f"{base_wcs_url}&REQUEST=ProcessCoverages&QUERY={query}"
                print("                        ")
                print("                        ")
                print(url)
                print("                        ")
                print("                        ")
                response = requests.get(url, auth=(rasdaman_username, rasdaman_password))
                #valls=(str(response.text)[1:-1])
                #sample_result.append(ansi_val)
                if response.status_code == 200:
                    sample_result.append(valls)
                    print("RESULT:", valls)
                    print("   ")
                    print("_______________________________________")
                else:
                    sample_result.append(response.status_code)
            else:
                print("NOT COVERED GEOGRAPHICALLY:")
                print("_______________________________________")
                #print("Xmax=", xmax, "YMAX=", ymax, "Xmin", xmin, "Ymin", ymin)
                #print("Xquery=", x1, "Yquery=", y1)
                #print(ansi_val, time_min, time_max)
                valls=("NOT COVERED")
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




def requestDataWGS(infoheader,layerlist,samples, filepath):
    env_vars = dotenv_values()
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?SERVICE=WCS&VERSION=2.0.1"
    result=[]
    header=[]
    header.append('sample,')
    header.append('lat,')
    header.append('lon,')
    header.append('date_slice,')
    if len(layerlist) > 1:
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry][0]
            header.append(layer)
    else:
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry]
            header.append(layer)
    log=open("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/WormpickerRequestWGS.log", "a")
    sys.stdout = log
    all_results=[] 
    for key, value in samples.items():
        sample_result=[]
        latitude_sample=value[0] 
        longitude_sample=value[1]
        name_date = key
        print(" ")
        print(" ")
        print("NEW SAMPLE")
        print("COORDINATES:", latitude_sample, longitude_sample)
        try:
            date = value[2]
            #print("DATE ASKED:",date)
            type(date)
        except IndexError:
            date = str(name_date[-10:])
            #print("DATE ASKED:",date)
        try:
            daydate=round_down_to_day(date)
            converted_date=datetime.strptime(daydate, "%Y-%m-%dT%H:%M:%S").isoformat() + "Z"
            time_asked=converted_date
            #print("correct format")
        except ValueError:
            date = process_date_string(date + 'T00:00:00Z')
            time_asked=str(date)
            #print(time_asked, name_date)
            #print(sample_result)
        sample_result.append(name_date)
        sample_result.append(latitude_sample)
        sample_result.append(longitude_sample)
        sample_result.append(time_asked)
        all_results.append(sample_result)
        #print(sample_result)
        print("QUERYING:", time_asked)
        try:
            converted_date = datetime.strptime(date, "%Y-%m-%d").isoformat() + "Z"
            search_time=parser.isoparse(converted_date)
            #print(search_time)
        except ValueError:
            #print("converted_cate=date")
            search_time=parser.isoparse(time_asked)
        for data_entry in range(0,len(layerlist)):
            i_header=infoheader
            #print(i_header)
            latitude_request = latitude_sample
            longitude_request = longitude_sample
            layer=layerlist[data_entry][i_header.index("CoverageID")]
            crs=layerlist[data_entry][i_header.index("CRS")]
            time_min= layerlist[data_entry][i_header.index("f")]
            time_max= layerlist[data_entry][i_header.index("t")]
            time_min=time_min.replace('"','')
            time_max=time_max.replace('"','')
            datetime1 = parser.isoparse(time_min)
            datetime2 = parser.isoparse(time_max)
            AxisLabels=layerlist[data_entry][i_header.index("axislabels")]
            TimeLabel=AxisLabels.split(",")[0]
            YLabel=AxisLabels.split(",")[1]
            XLabel=AxisLabels.split(",")[2]
            #print(TimeLabel,YLabel, XLabel)
            ymin= float(layerlist[data_entry][i_header.index("minlat")])
            ymax = float(layerlist[data_entry][i_header.index("maxlat")])
            xmin= float(layerlist[data_entry][i_header.index("minlong")])
            xmax = float(layerlist[data_entry][i_header.index("maxlong")])
            xres=float(layerlist[data_entry][i_header.index("resolution_dim1")])
            yres=float(layerlist[data_entry][i_header.index("resolution_dim2")])
            #strlat="&subset=Y({},{})"
            #strlon = "&subset=X({},{})"
            ansi_str_o="&subset=ansi(\"{}\")"
            ansi_str_n=ansi_str_o.replace("ansi",TimeLabel)
            if search_time > datetime1 and search_time < datetime2:
                #print(search_time, "- COVERED TEMPROALLY")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                ##produce timestamp here
                ansi_val=time_asked
            else:
                min_date = datetime1 if abs((datetime1 - search_time).days) < abs((datetime2 - search_time).days) else datetime2
                dist=abs(min_date-search_time).days
                #print(time_asked, "NOT COVERED TEMPORALLY, querying:", min_date, ". This is", dist, "days apart!")
                #set timestamp here a general one, the one were coverage is there
                ansi_val=time_min
            #print("CRS of layer:", crs)
            #print("Calculating Boudns of Layer for EPSG:4326:")
            y1=float(latitude_request)
            x1=float(longitude_request)
            #strlat="&subset=Lat({},{})"
            #strlon = "&subset=Long({},{})"
            #print(layer)
            crs_indicator= crs.replace("EPSG/0/", "EPSG:")
            output_format = "text/csv"
            #print(name_date, request_cov_id, request_encode_format, layer)
            #min_y, min_x = trans4mEPSG("EPSG:3035","EPSG:4326", float(ymin), float(xmin))
            #max_y, max_x = trans4mEPSG("EPSG:3035","EPSG:4326",float(ymax), float(xmax))
            #print("Geobounds LAT:", ymin, "-", ymax)
            #print("Geobounds LONG:", xmin, "-", xmax)
            if x1 <= xmax and x1 >= xmin and y1 >= ymin and y1 <= ymax and ansi_val != time_min and ansi_val != time_max:
                #print("The Layer:", layer, " is within temporal coverage. Querying for", time_asked, "as", ansi_val)
                #TimeLabel="date"
                #XLabel="Lon"
                #YLabel="Lat"
                query = f"for $c in ({layer}) return encode($c[{TimeLabel}(\"{ansi_val}\"),{XLabel} :\"EPSG:4326\"( {longitude_sample} ), {YLabel}:\"EPSG:4326\"( {latitude_sample}) ], \"{output_format}\")"
                #query = f"for $c in ({layer}) return encode($c[date:(\"{ansi_val}\"),Lon :\"CRS:1\"{grid_indices_axis_X},Lat :\"CRS:1\"{grid_indices_axis_Y} ], \"{output_format}\")"
                # Construct the URL with variables
                url = f"{base_wcs_url}&REQUEST=ProcessCoverages&QUERY={query}"
                #print("                        ")
                #print("                        ")
                print(url)
                #print("                        ")
                print("                        ")
                response = requests.get(url, auth=(rasdaman_username, rasdaman_password))
                valls=str(response.text)
                #sample_result.append(ansi_val)
                if response.status_code == 200:
                    sample_result.append(valls)
                    print("RESULT:", valls)
                    print("   ")
                    print("_______________________________________")
                else:
                    sample_result.append(response.status_code)
            else:
                #print("NOT COVERED GEOGRAPHICALLY:")
                #print("_______________________________________")
                #print("Xmax=", xmax, "YMAX=", ymax, "Xmin", xmin, "Ymin", ymin)
                #print("Xquery=", x1, "Yquery=", y1)
                #print(ansi_val, time_min, time_max)
                valls=("NOT COVERED")
                sample_result.append(valls)
        #sample_result.append(query)        
        result.append(sample_result)
        #header.append('query_sent,')
        
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