import sys
import os 
import csv
import requests
from dotenv import dotenv_values
from dateutil import parser     
from datetime import datetime, timedelta
#from module_crs_converter import trans4mEPSG
import math
from osgeo import gdal, osr
import xmltodict
import re

import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox

#env_path="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/.env"

## Load environment variables from the file
#with open(env_path) as f:
#    for line in f:
#        key, value = line.strip().split('=', 1)
#        os.environ[key] = value

# Access logs later


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


#####
def remove_partial_dates(data_dict):
    # Regular expression pattern for "MM-DDTT00"
    invalid_data_pattern = re.compile(r"\d{4}-\d{2}-\d{2}$")
    # List to store keys to remove
    keys_to_remove = []
    # # Identify keys with partial dates
    for key in data_dict.keys():
        date_part= key.split('_')[-1]
        if not invalid_data_pattern.match(date_part):
            keys_to_remove.append(key)
            # Remove identified keys 
    for key in keys_to_remove:
        del data_dict[key]
    return data_dict


#def is_valid_iso_date(date_string):
#    try:
#        # Attempt to parse the string as an ISO-formatted date
#        datetime.fromisoformat(date_string)
#        return True
#    except ValueError:
#        return False
#
#def process_date_string(date_string):
#    if is_valid_iso_date(date_string):
#        print(f"{date_string} is already in ISO format. Doing nothing.")
#        return date_string
#    else:
#        # Add your processing logic here if the date_string is not in ISO format
#        print(f"Processing {date_string}...")
#        # For now, just returning the input date_string as an example
#        return date_string
#
#def get_grid_indices_for_axis_X(geo_lower_bound, geo_upper_bound, xmin, xmax, xres):       
#    print("geo_lower_bound_X: {}, geo_upper_bound_X: {}".format(geo_lower_bound, geo_upper_bound))
#    lower_grid_index = math.floor( (geo_lower_bound - xmin) / xres )
#    upper_grid_index = math.floor( (geo_upper_bound - xmin) / xres )
#    if geo_upper_bound == xmax:
#        upper_grid_index = upper_grid_index - 1       
#    return "({}:{})".format(lower_grid_index, upper_grid_index)
#    
#def get_grid_indices_for_axis_Y(geo_lower_bound, geo_upper_bound, ymin, ymax, yres):
#    print("geo_lower_bound_Y: {}, geo_upper_bound_Y: {}".format(geo_lower_bound, geo_upper_bound))
#    lower_grid_index = math.floor( (geo_upper_bound - ymax) / yres )
#    upper_grid_index = math.floor( (geo_lower_bound - ymax) / yres )
#    if geo_lower_bound == ymin:
#        upper_grid_index = upper_grid_index - 1   
#    return "({}:{})".format(lower_grid_index, upper_grid_index)
#
#
#def round_down_to_day(date_str):
#    # Parse the date string
#    try:
#        date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S")
#    except ValueError:
#        try:
#            date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H:%M")
#        except ValueError:
#            date_obj = datetime.strptime(date_str, "%Y-%m-%dT%H")
#    rounded_down_date = date_obj.replace(hour=0, minute=0, second=0)
#    rounded_down_date_str = rounded_down_date.strftime("%Y-%m-%dT%H:%M:%S")
#    return rounded_down_date_str


def findDate(inputtime, description):
    # Example format of your datetime strings, adjust as necessary
    datetime_format = "%Y-%m-%dT%H:%M:%S.%fZ"
    date_obj = datetime.strptime(inputtime, "%Y-%m-%d")
    timestring = date_obj.strftime("%Y-%m-%dT%H:%M:%S.%fZ")[:-4] + "Z"
    testtime = datetime.strptime(timestring, datetime_format)
    min_diff_days = float('inf')
    withinTime=False
    closest_date = None
    try:
        b=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:metadata']['gmlcov:Extension']['ras:covMetadata']['ras:axes']['ras:time']['ras:areasOfValidity']['ras:area']
        for item in b: 
            starttime=datetime.strptime(item['@start'], datetime_format)
            endtime=datetime.strptime(item['@end'], datetime_format)
            if (starttime < testtime < endtime):
                withinTime=True
                closest_date=testtime
                min_diff_days=0
                return closest_date,min_diff_days
                break
            else:
                start_diff_days = abs((starttime - testtime).days)
                end_diff_days = abs((endtime - testtime).days)
                if start_diff_days < min_diff_days:
                    min_diff_days = start_diff_days
                    closest_date = starttime + timedelta(days=1)
                if end_diff_days < min_diff_days:
                    min_diff_days = end_diff_days
                    closest_date = endtime - timedelta(days=1)
                if not withinTime and closest_date is not None:
                    continue
                    #print("Sample time not within boundaries. Finding closest time slice....")
        return closest_date, min_diff_days
    except KeyError:
        #print("Trying to find right key.")
        starttime_str=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['gml:lowerCorner'].split(" ")[0][1:-1]
        endtime_str=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['gml:upperCorner'].split(" ")[0][1:-1]
        starttime=datetime.strptime(starttime_str, datetime_format)
        endtime=datetime.strptime(endtime_str, datetime_format)
        if starttime < testtime < endtime:
            #print("Sample time", testtime, "within time ranges of layer:", starttime_str, endtime_str)
            return testtime, "0"
        else:
            start_diff_days = abs((starttime - testtime).days)
            end_diff_days = abs((endtime - testtime).days)
            if start_diff_days < min_diff_days:
                min_diff_days = start_diff_days
                closest_date = starttime + timedelta(days=1)
            if end_diff_days < min_diff_days:
                min_diff_days = end_diff_days
                closest_date = endtime - timedelta(days=1)
            return closest_date,min_diff_days
    except TypeError:
        #print("Handling Type Error.")
        try:
            timestamps=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gmlrgrid:ReferenceableGridByVectors']['gmlrgrid:generalGridAxis'][0]['gmlrgrid:GeneralGridAxis']['gmlrgrid:coefficients']
            if timestring in timestamps:
                #print(timestring, "Sample time within time ranges of layer.")
                return timestring, "0"
            else:
                for time in timestamps.split(" "):
                    #print(time)
                    #print(testtime)
                    diff = abs((testtime - datetime.strptime(time[1:-1], datetime_format)).days)
                    if diff < min_diff_days:
                        min_diff_days = diff
                        closest_date = time[1:-1]
                return closest_date, min_diff_days
        except KeyError:
            starttime_str=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['gml:lowerCorner'].split(" ")[0][1:-1]
            endtime_str=description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['gml:upperCorner'].split(" ")[0][1:-1]
            #print(starttime_str, endtime_str)
            starttime=datetime.strptime(starttime_str, datetime_format)
            endtime=datetime.strptime(endtime_str, datetime_format)
            if starttime < testtime < endtime:
                #print("Sample time within time ranges of layer.",)
                return testtime, "0"
            else:
                start_diff_days = abs((starttime - testtime).days)
                end_diff_days = abs((endtime - testtime).days)
                if start_diff_days < min_diff_days:
                    min_diff_days = start_diff_days
                    closest_date = starttime + timedelta(days=1)
                if end_diff_days < min_diff_days:
                    min_diff_days = end_diff_days
                    closest_date = endtime - timedelta(days=1)
                return closest_date,min_diff_days
###wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['gml:upperCorner'].split(" ")[0][1:-1]
#result,diff=findDate("2021-03-01",wcs_coverage_description)



def requestDataWGS(infoheader,layerlist,samples, filepath,offset=0, approximate=True, rasdaman_username=None, rasdaman_password=None, rasdaman_endpoint=None, loggerobject=None):
    #sys.stdout = log
    datetime.now()
    env_vars = dotenv_values()
    #rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT') 
    #rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    #rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    rasdaman_endpoint = rasdaman_endpoint or env_vars.get('RASDAMAN_SERVICE_ENDPOINT') 
    rasdaman_username = rasdaman_username or env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = rasdaman_password or env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?SERVICE=WCS&VERSION=2.1.0"
    result=[]
    distances=[]
    header=[]
    header.append('sample')
    header.append('lat')
    header.append('lon')
    for data_entry in range(0,len(layerlist)):
            i_header=infoheader
            if len(layerlist) > 1:
                number_bands=len(layerlist[data_entry][-2])
                layer=layerlist[data_entry][0]
                layerbands_nr=layerlist[data_entry][-2]
                for offday in range(abs(offset)+1):
                    for col in layerbands_nr:
                        if offset > 0:
                            colname_new=str(layer)+"_"+str(col)+"+"+str(offday)
                        if offset < 0:
                            colname_new=str(layer)+"_"+str(col)+"-"+str(offday)
                        else:
                            colname_new=str(layer)+"_"+str(col)
                        header.append(colname_new)
    for key, value in samples.items():
        sample_result=[]
        sample_distances=[]
        latitude_sample=value[0] 
        longitude_sample=value[1]
        name_date = key
        loggerobject.info("PROCESSING NEW SAMPLE")
        loggerobject.info(f"COORDINATES - Input Latitude: {latitude_sample}, Input Longitude {longitude_sample}")
        try:
            date = value[2]
            loggerobject.info(f"DATE SAMPLE: {date}")
            type(date)
        except IndexError:
            #date = str(name_date[-10:])
            #loggerobject.warning(f"DATE QUERIED taken from samplename: {date} ")
            loggerobject.warning("No date for sample given. Please provide a date column in the input file.")
            return
        sample_result.append(name_date)
        sample_result.append(latitude_sample)
        sample_result.append(longitude_sample)
        #print("TRYING:", date)
        #days_offset=0+offset
        for data_entry in range(0,len(layerlist)):
            value_result=[]
            i_header=infoheader
            ##old
            #print(i_header)
            if len(layerlist) > 1:
                number_bands=len(layerlist[data_entry][-2])
                #print(number_bands)
                ##new block 2 lines
                layer=layerlist[data_entry][0]
                layerbands_nr=layerlist[data_entry][-2]
                for offday in range(abs(offset)+1):
                    for col in layerbands_nr:
                        if offset > 0:
                            colname_new=str(layer)+"_"+str(col)+"+"+str(offday)
                            #sample_result.append([colname_new])
                        if offset < 0:
                            colname_new=str(layer)+"_"+str(col)+"-"+str(offday)
                        else:
                            colname_new=str(layer)+"_"+str(col)
                            #sample_result.append([colname_new])
                        #header.append(colname_new)
                        value_result.append(colname_new)
            else:
                number_bands=len(layerlist[0][-2])
                layer=layerlist[0][0]
                layerbands_nr=layerlist[0][-2]
                for offday in range(abs(offset)+1):
                    for col in layerbands_nr:
                        if offset >= 0:
                            colname_new=str(layer)+"_"+str(col)+"+"+str(offday)
                            #sample_result.append([colname_new])
                        else:
                            colname_new=str(layer)+"_"+str(col)+"-"+str(offday)
                            #sample_result.append([colname_new])
                        value_result.append(colname_new)
            layer=layerlist[data_entry][i_header.index("CoverageID")]
            loggerobject.info(f"")
            loggerobject.info(f"Processing Coverage {layer}")
            ###add description URL and get time stamps
            describe_url='https://fairicube.rasdaman.com/rasdaman/ows?&SERVICE=WCS&VERSION=2.1.0'
            response = requests.get(describe_url + "&REQUEST=DescribeCoverage&COVERAGEID=" + layer,auth=(rasdaman_username, rasdaman_password))
            try:
                if 'wcs:CoverageDescriptions' in xmltodict.parse(response.content).keys():
                    wcs_coverage_description = xmltodict.parse(response.content)
                elif 'ows:ExceptionReport' in xmltodict.parse(response.content).keys():
                    loggerobject.debug("Credentials provided are not valid.")
                #print(response.content)
                #print(wcs_coverage_description)
            except:
                loggerobject.debug(f"Response Content: {response.content}")
                continue
            approxvar=approximate
            if approxvar is True:
                loggerobject.info(f"Trying to find the date {date}")
                try:
                    date_obj = datetime.strptime(date, "%Y-%m-%d")
                except:
                    return
                try:
                    timetoquery,differencedays=findDate(date,wcs_coverage_description)
                    #date_obj = datetime.strptime(timetoquery_raw, "%Y-%m-%d")
                    #timetoquery = timetoquery_raw.strftime("%Y-%m-%dT%H:%M:%S.%fZ")[:-4] + "Z"
                    #print(timetoquery)
                    loggerobject.info(f"Querying approximate time: {timetoquery}")
                    loggerobject.info(f"Difference in Days: {differencedays}")
                except:
                    loggerobject.error("Could not execute function called: findDate. Either date or coverade_description cannot be assigned." )
                    sample_result.append("NA")
                    continue
            else:
                date_obj = datetime.strptime(date, "%Y-%m-%d")
                timetoquery = date_obj.strftime("%Y-%m-%dT%H:%M:%S.%fZ")[:-4] + "Z"
                differencedays="0"
                loggerobject.warning(f"Querying sample date exactly without approximation at {timetoquery} . This can lead to sparse data. ")
            AxisLabels=layerlist[data_entry][i_header.index("axislabels")]
            if AxisLabels == "ds.earthserver.xyz":
                loggerobject.info("Querying Earth Server Data")
                #TimeLabel="ansi"
                TimeLabel=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['@axisLabels'].split(" ")[0]
                #YLabel="lat"
                YLabel=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['@axisLabels'].split(" ")[1]
                #XLabel="lon"
                XLabel=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:boundedBy']['gml:Envelope']['@axisLabels'].split(" ")[2]
            else:
                TimeLabel=AxisLabels.split(",")[0]
                YLabel=AxisLabels.split(",")[1]
                XLabel=AxisLabels.split(",")[2]
            loggerobject.info(f"Axislabels: {TimeLabel} {YLabel} {XLabel}")
            #ansi_str_o="&subset=ansi(\"{}\")"
            ###### NEW TIME LABEL
            output_format = "text/csv"
            for dayoff in range(abs(offset)+1):
                #date_obj = datetime.strptime(timetoquery, "%Y-%m-%dT%H:%M:%S.%fZ")
                if offset >= 0 and isinstance(timetoquery, str):
                    loggerobject.info("Approximation parameter set. Offset is positive.")
                    date2=datetime.strptime(timetoquery[0:10], "%Y-%m-%d") + timedelta(days=dayoff)
                    #print(date2)
                elif offset < 0 and isinstance(timetoquery, str):
                    date2=datetime.strptime(timetoquery[0:10], "%Y-%m-%d") - timedelta(days=dayoff)
                elif offset >= 0 and isinstance(timetoquery, datetime):
                    date2=timetoquery + timedelta(days=dayoff)
                else:
                    date2=timetoquery - timedelta(days=dayoff)
                timetoquery = date2.strftime("%Y-%m-%dT%H:%M:%S.%fZ")[:-4] + "Z"
                query = f"for $c in ({layer}) return encode($c[{TimeLabel}(\"{timetoquery}\"),{XLabel} :\"EPSG:4326\"( {longitude_sample} ), {YLabel}:\"EPSG:4326\"( {latitude_sample}) ], \"{output_format}\")"
                url = f"{base_wcs_url}&REQUEST=ProcessCoverages&QUERY={query}"
                loggerobject.info(f" START QUERYING: {url}")
                response = requests.get(url, auth=(rasdaman_username, rasdaman_password))
                value_response=str(response.text)
                #print(value_response)
                if isinstance(value_response, str) and value_response.startswith("{") and value_response.endswith("}"):
                    valls=str(response.text)[1:-1].replace(" ", ",")
                else:
                    valls=str(response.text).replace(" ", ",")
                if response.status_code == 200 and valls != '':
                    #sample_result.append(valls)
                    for singleval in valls.strip('"').split(","):
                        #print(singleval)
                        #singlevalq="{"+singleval+"}"
                        sample_result.append(singleval)
                        sample_distances.append(differencedays)
                        value_result.append(singleval)
                    #sample_result.append(layer)
                    loggerobject.info(f"RESULT: {valls}")
                    loggerobject.info("_______________________________________")
                else:
                    loggerobject.info(f" {response}")
                    for _ in range(number_bands):
                        #t="NA"+str(_)
                        t="NA"
                        sample_result.append(t)
                        sample_distances.append("NA")
                        value_result.append("NA")
                loggerobject.info(f" {value_result}")
                    #sample_result.append(layer)
        result.append(sample_result)
        distances.append(sample_distances)
        #result.append(samples_results_sep)
        #header.append('query_sent,')  
    #sys.stdout = sys.__stdout__
    datetime.now()
    sys.stdout = sys.__stdout__
    if filepath!="NONE":
        loggerobject.info(f"RESULT WRITTEN TO OUTPUT FILE: {filepath}")
        with open(filepath, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            writer.writerows(result)
        #with open((filepath+"2"), "w", newline="") as csvfile:
        #    writer = csv.writer(csvfile)
        #    writer.writerow(header[3:])
        #    writer.writerows(distances)
            #writer.writerows(value_result)
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