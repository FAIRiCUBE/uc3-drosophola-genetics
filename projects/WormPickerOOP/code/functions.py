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

#env_path="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/.env"

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


import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox


mode="manual"


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





def requestDataWGS(infoheader,layerlist,samples, filepath, approximate="FALSE"):
    env_vars = dotenv_values()
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?SERVICE=WCS&VERSION=2.1.0"
    result=[]
    header=[]
    header.append('sample,')
    header.append('lat,')
    header.append('lon,')
    if len(layerlist) > 1:
        for data_entry in range(0,len(layerlist)):
            layer=layerlist[data_entry][0]
            layerbands_nr=layerlist[data_entry][-1]
            #header.append(layer)
            for col in layerbands_nr:
                colname_new=str(layer)+str(col)
                header.append(colname_new)
    else:
        #for data_entry in range(0,len(layerlist)):
        layer=layerlist[0][0]
        layerbands_nr=layerlist[0][-1]
        print(layerbands_nr)
        for col in layerbands_nr:
            colname_new=str(layer)+"_"+str(col)
            header.append(colname_new)
            #header.append
            #header.append(layer)
    log=open("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/LatestResults.log", "a")
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
        print("COORDINATES lat/long:", latitude_sample, longitude_sample)
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
        #all_results.append(sample_result)
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
            if len(layerlist) > 1:
                number_bands=len(layerlist[data_entry][-1])
            else:
                number_bands=len(layerlist[0][-1])
            latitude_request = latitude_sample
            longitude_request = longitude_sample
            layer=layerlist[data_entry][i_header.index("CoverageID")]
            print(layer)
            crs=layerlist[data_entry][i_header.index("CRS")]
            time_min= layerlist[data_entry][i_header.index("f")]
            time_max= layerlist[data_entry][i_header.index("t")]
            time_min=time_min.replace('"','')
            time_max=time_max.replace('"','')
            datetime1 = parser.isoparse(time_min)
            datetime2 = parser.isoparse(time_max)
            AxisLabels=layerlist[data_entry][i_header.index("axislabels")]
            if AxisLabels == "ds.earthserver.xyz":
                TimeLabel="time"
                YLabel="lat"
                XLabel="lon"
            else:
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
                #print(search_time, "- COVERED TEMPROALLY")  #make this to acommand which produces a file, saying which samples and layers REALLYY overlap (also in time)
                ##produce timestamp here
                ansi_val=time_asked
            elif search_time < datetime1  and approximate=="TRUE" or search_time > datetime2 and approximate=="TRUE":
                min_date = datetime1 if abs((datetime1 - search_time).days) < abs((datetime2 - search_time).days) else datetime2
                dist=abs(min_date-search_time).days
                print(time_asked, "NOT COVERED TEMPORALLY, querying:", min_date, ". This is", dist, "days apart!")
                #set timestamp here a general one, the one were coverage is there
                ansi_val=min_date
            else:
                print("Something is wrong with temporal approximation!")
            #print("CRS of layer:", crs)
            #print("Calculating Boudns of Layer for EPSG:4326:")
            y1=float(latitude_request)
            x1=float(longitude_request)
            #strlat="&subset=Lat({},{})"
            #strlon = "&subset=Long({},{})"
            #print(layer)
            crs_indicator= crs.replace("EPSG/0/", "EPSG:")
            output_format = "text/csv"
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
                value_response=str(response.text)
                print(value_response)
                if isinstance(value_response, str) and value_response.startswith("{") and value_response.endswith("}"):
                    valls=str(response.text)[1:-1].replace(" ", ",")
                else:
                    valls=str(response.text).replace(" ", ",")
                #if isinstance(value_response, tuple):
                    #valls=value_response[0]
                #sample_result.append(ansi_val)
                if response.status_code == 200:
                    #sample_result.append(valls)
                    for singleval in valls.strip('"').split(","):
                        print(singleval)
                        sample_result.append(singleval)
                    print("RESULT:", valls)
                    print("   ")
                    print("_______________________________________")
                else:
                    #sample_result.append(response.status_code)
                    #sample_result.append("na")
                    #na_adjust= ["na"] * number_bands
                    #sample_result.append(na_adjust)
                    for _ in range(number_bands):
                        sample_result.append("na")
            else:
                #print("NOT COVERED GEOGRAPHICALLY:")
                #print("_______________________________________")
                #valls=["nc"] * number_bands
                #sample_result.append(valls)
                for _ in range(number_bands):
                        sample_result.append("nc")
        #sample_result.append(query)   
        result.append(sample_result)
        #result.append(samples_results_sep)
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
    



####### 

def summarize_csv(file_path,output_text_file):
    # Initialize a dictionary to hold the summary
    summary = {}
    column_data = {}
    # Open and read the CSV file
    with open(file_path, mode='r') as csvfile:
        reader = csv.reader(csvfile)      
        # Read header
        headers = next(reader)        
        # Initialize columns in the summary dictionary
        for header in headers:
            summary[header] = {
                'NA Count': 0,
                'NC Count': 0,
                'Non-Zero Numeric Count': 0,
                'Total Count': 0
            }
            column_data[header] = []       
        # Process each row
        for row in reader:
            for i, value in enumerate(row):
                header = headers[i]
                column_data[header].append(value)
                summary[header]['Total Count'] += 1 
                if value == "na":
                    summary[header]['NA Count'] += 1
                elif value == "nc":
                    summary[header]['NC Count'] += 1
                else:
                    try:
                        # Convert to float and check if it is not zero
                        numeric_value = float(value)
                        if numeric_value != 0:
                            summary[header]['Non-Zero Numeric Count'] += 1
                    except ValueError:
                        # Ignore values that cannot be converted to float
                        pass
    # Prepare to write summary to text file
    with open(output_text_file, 'w') as file:
        for column, stats in summary.items():
            total_count = stats['Total Count']
            na_percentage = (stats['NA Count'] / total_count) * 100 if total_count > 0 else 0
            nc_percentage = (stats['NC Count'] / total_count) * 100 if total_count > 0 else 0
            non_zero_numeric_percentage = (stats['Non-Zero Numeric Count'] / total_count) * 100 if total_count > 0 else 0
            
            # Write summary for each column
            file.write(f"Column: {column}\n")
            file.write(f"  NA Percentage: {na_percentage:.2f}%\n")
            file.write(f"  NC Percentage: {nc_percentage:.2f}%\n")
            file.write(f"  Non-Zero Numeric Percentage: {non_zero_numeric_percentage:.2f}%\n")
            file.write("\n")