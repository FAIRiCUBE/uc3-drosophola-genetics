import tkinter as tk
import tkinter as ttk
from tkinter import simpledialog
import requests
import xmltodict
import csv
from pathlib import Path
from dotenv import dotenv_values
import re

def getLayers(savepath="NONE"):
    #create an empty array and add headers
    layer_info_2 = []
    layer_info_2.append(("CoverageID", "CRS", "minlong","maxlong", "minlat", "maxlat","f","t", "axislabels", "resolution_dim1","resolution_dim2" , "resolution_time" , "bands", "null_values"))
    #Load environment variables from the .env file
    env_vars = dotenv_values()
    # Access specific environment variables (useful when declared in a unix environment already)
    rasdaman_endpoint = env_vars.get('RASDAMAN_SERVICE_ENDPOINT')
    rasdaman_username = env_vars.get('RASDAMAN_CRED_USERNAME')
    rasdaman_password = env_vars.get('RASDAMAN_CRED_PASSWORD')
    base_wcs_url = rasdaman_endpoint + "?&SERVICE=WCS&ACCEPTVERSIONS=2.1.0"
    response = requests.get(base_wcs_url + "&REQUEST=GetCapabilities", auth=(rasdaman_username, rasdaman_password))
    wcs_capabilities = xmltodict.parse(response.content)
    wcs_coverage_summary = wcs_capabilities['wcs20:Capabilities']['wcs20:Contents']['wcs20:CoverageSummary']
    #print(json.dumps(wcs_coverage_summary, indent=2))
    type(wcs_coverage_summary)
    #cov_id_list=[]
    #extract all coverage_ids and informations from the service endpoint (working on getCapabilities results)
    describe_url='https://fairicube.rasdaman.com/rasdaman/ows?&SERVICE=WCS&VERSION=2.1.0'
    for i in range(0,len(wcs_coverage_summary)):
        coverage_id= wcs_coverage_summary[i]['wcs20:CoverageId']
        coverageSubtype=wcs_coverage_summary[i]['wcs20:CoverageSubtype']
        dim=int(wcs_coverage_summary[i]['ows:BoundingBox']['@dimensions'])
        #crs_str = wcs_coverage_summary[i]['ows:BoundingBox']['@crs'][-11:]
        crs_str = wcs_coverage_summary[i]['ows:BoundingBox']['@crs']
        crs = '/'.join(crs_str.split('/')[-3:])
        if dim==3:
            print("DIM IS 3", coverage_id)
            x_min=wcs_coverage_summary[i]['ows:WGS84BoundingBox']['ows:LowerCorner'].split(" ")[0]
            y_min=wcs_coverage_summary[i]['ows:WGS84BoundingBox']['ows:LowerCorner'].split(" ")[1]
            x_max=wcs_coverage_summary[i]['ows:WGS84BoundingBox']['ows:UpperCorner'].split(" ")[0]
            y_max=wcs_coverage_summary[i]['ows:WGS84BoundingBox']['ows:UpperCorner'].split(" ")[1]
            bb_low=wcs_coverage_summary[i]['ows:BoundingBox']['ows:LowerCorner']
            bb_upp=wcs_coverage_summary[i]['ows:BoundingBox']['ows:UpperCorner']
            #x_min=bb_low.split(" ")[1] 
            #y_min=bb_low.split(" ")[2]
            date_min=bb_low.strip().split(" ")[0][1:-1]
            date_max=bb_upp.strip().split(" ")[0][1:-1]
            #x_max=bb_upp.split(" ")[1] 
            #y_max=bb_upp.split(" ")[2]
            try:
                axislabels=wcs_coverage_summary[i]['ows:AdditionalParameters']['ows:AdditionalParameter'][3]['ows:Value']
                if re.search("ds.earthserver.xyz", axislabels):
                    print("matched: ", axislabels)
                    axislabels=wcs_coverage_summary[i]['ows:AdditionalParameters']['ows:AdditionalParameter'][1]['ows:Value']
                    if re.search(r"[0-9]", axislabels):
                        axislabels=wcs_coverage_summary[i]['ows:AdditionalParameters']['ows:AdditionalParameter'][2]['ows:Value']
                        print("substituted with:", axislabels)
            except IndexError:
                try:
                    axislabels=wcs_coverage_summary[i]['ows:AdditionalParameters']['ows:AdditionalParameter'][2]['ows:Value']
                except IndexError:
                    axislabels=wcs_coverage_summary[i]['ows:AdditionalParameters']['ows:AdditionalParameter'][1]['ows:Value']
            #print(coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_min)
            #!!!!usually:layer_info_2.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_max,axislabels))
            #try:
                #for ID in range(1,len(layer_info_2)):
                    #coverage=layer_info_2[ID][0]
            print(coverage_id)
            response = requests.get(describe_url + "&REQUEST=DescribeCoverage&COVERAGEID=" + coverage_id,auth=(rasdaman_username, rasdaman_password))
            wcs_coverage_description = xmltodict.parse(response.content)
            null_values=[]
            try:
                try:
                    null_value=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field']['swe:Category']['swe:nilValues']['swe:NilValues']['swe:nilValue']['#text']
                except KeyError:
                    try:
                        null_value=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field']['swe:Quantity']['swe:nilValues']['swe:NilValues']['swe:nilValue']['#text']
                    except TypeError:
                        null_value=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field']['swe:Quantity']['swe:nilValues']['swe:NilValues']['swe:nilValue'][0]['#text']
            except:
                null_value="NA"       
        #print(json.dumps(wcs_coverage_description, indent=2))
            null_values.append(null_value)
            rr=[]
            try: 
                l=len(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gmlrgrid:ReferenceableGridByVectors']['gmlrgrid:generalGridAxis'][0]['gmlrgrid:GeneralGridAxis']['gmlrgrid:offsetVector']['#text'].split(" "))
                for i in reversed(range(0,l)):
                    resi=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gmlrgrid:ReferenceableGridByVectors']['gmlrgrid:generalGridAxis'][i]['gmlrgrid:GeneralGridAxis']['gmlrgrid:offsetVector']['#text'] 
                    resolution=resi.split(" ")[i]
                    rr.append(resolution)
                    #print(resolution)
            except KeyError:
                xml_info=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gml:domainSet']['gml:RectifiedGrid']['gml:offsetVector']
                l=len(xml_info)
                for i in range(0,l):
                    resolution=xml_info[i]['#text'].split(" ")[i]
                    rr.append(resolution)
            band_infos=[]            
            try:
                if isinstance(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'],list):
                    bands=len(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'])
                    for band_nr in range(0,bands): 
                        band_name=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'][band_nr]['@name']
                        print(band_name)
                        band_infos.append(band_name)
                        try:
                            null_value=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'][band_nr]['swe:Quantity']['swe:nilValues']['swe:NilValues']['swe:nilValue']['#text']
                        except TypeError:
                            null_value=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'][band_nr]['swe:Quantity']['swe:nilValues']['swe:NilValues']['swe:nilValue']
                        null_values.append(null_value)
                    #print(band_infos)
                if isinstance(wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field'],dict):
                    band_name=wcs_coverage_description['wcs:CoverageDescriptions']['wcs:CoverageDescription']['gmlcov:rangeType']['swe:DataRecord']['swe:field']['@name']
                    band_infos.append(band_name)              
            except KeyError:
                pass
            #x=tuple([resolution])
        #print(layer_info_2[ID])
            layer_info_2.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_max,axislabels, rr[0], rr[1], rr[2], band_infos, null_values))
        else:
            print("NOT ALL DIMENSIONS AVAILABLE")
        #layer_info_2.append((coverage_id, crs, x_min, x_max, y_min, y_max, date_min, date_max,axislabels, rr[0], rr[1], rr[2], band_infos))
    if savepath!="NONE":
        os.chdir(savepath)
        with open(savepath, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(layer_info_2)
    else:
        return layer_info_2

def processFunctionResult(result):
    global coverage_instance, boundary_button
    coverage_instance = Coverage(result)
    print("Coverage instance created.")
    boundary_button.config(state=tk.NORMAL)  # Enable the Boundary button

def processCoverageInstance(result):
    global boundary
    global samplescovered
    boundary = result.getBoundary()
    result.getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")
    samplescovered = result.samples
    print("Boundary created:", boundary, "Samples covered are:", samplescovered)


def execute_function():
    global function_result, proceed_button, boundary_button
    if savepath_var.get() == "NONE":
        function_result = getLayers(savepath="NONE")
        print("Function executed without custom savepath.")
    else:
        custom_savepath = simpledialog.askstring("Custom Savepath", "Enter custom savepath:")
        if custom_savepath:
            function_result = getLayers(savepath=custom_savepath)
            print("Function executed with custom savepath:", custom_savepath)
    proceed_button.config(state=tk.NORMAL)  # Enable the Proceed button
    proceed_button.focus_set()  # Set focus to the Proceed button
    execute_button.config(state=tk.DISABLED)  # Disable the Execute button

def proceed_function():
    global function_result, coverage_instance, proceed_button, boundary_button
    if 'function_result' in globals():
        processFunctionResult(function_result)
        proceed_button.config(state=tk.DISABLED)  # Disable the Proceed button
        execute_button.config(state=tk.DISABLED)  # Disable the Execute button
        boundary_button.config(state=tk.NORMAL)  # Enable the Boundary button
    else:
        print("Function result not available.")

def boundary_function():
    global coverage_instance, boundary_button
    if 'coverage_instance' in globals():
        processCoverageInstance(coverage_instance)
        boundary_button.config(state=tk.DISABLED)  # Disable the Boundary button
    else:
        print("Coverage instance not available.")

## Create the main application window
#root = tk.Tk()
#root.title("Function Executor")
#
## Create a label
#label = tk.Label(root, text="Select an option:")
#label.pack(padx=20, pady=10)
#
## Create a radio button for "NONE" option
#savepath_var = tk.StringVar(value="NONE")
#none_radio = tk.Radiobutton(root, text="Process with savepath=NONE", variable=savepath_var, value="NONE")
#none_radio.pack(anchor="w", padx=20)
#
## Create a radio button for custom savepath option
#custom_radio = tk.Radiobutton(root, text="Enter custom savepath", variable=savepath_var, value="CUSTOM")
#custom_radio.pack(anchor="w", padx=20)
#
## Create a button that executes the function
#execute_button = tk.Button(root, text="Execute Function", command=execute_function)
#execute_button.pack(padx=20, pady=10)
#
## Create a Proceed button (initially disabled)
#proceed_button = tk.Button(root, text="Make Coverage Object", state=tk.DISABLED, command=proceed_function)
#proceed_button.pack(padx=20, pady=10)
#
## Create a Boundary button (initially disabled)
#boundary_button = tk.Button(root, text="Get Coverage Boundary", state=tk.DISABLED, command=boundary_function)
#boundary_button.pack(padx=20, pady=10)
#
## Start the main event loop
#root.mainloop()
#
#layers_to_analyze=select_objects("manual",x)
#
#
#def select_layers_and_request_data_window():
#    x = coverage_instance
#    x.getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")
#    samplescovered = x.samples
#    layers_to_analyze = select_objects("manual", x)
#    layerlist = []
#    for layer in layers_to_analyze:
#        for i in range(0, len(x._data)):
#            entries_list = list(x._data[i])
#            try:
#                index = entries_list.index(layer)
#                layerlist.append(entries_list)
#            except:
#                continue
#    output_path = "/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/longtimetry_window.csv"
#    requestData(layerlist, samplescovered, output_path)
#    print("Processed layers and requested data.")
#
## Create the main application window
#root = tk.Tk()
#root.title("Select Layers and Request Data")
#
## Create a button for the process
#process_button = tk.Button(root, text="Process Layers and Request Data", command=select_layers_and_request_data_window)
#process_button.pack(padx=20, pady=10)
#
## Start the main event loop
#root.mainloop()