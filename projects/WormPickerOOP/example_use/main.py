from code.module_crs_converter import trans4mEPSG
from code.Objects import Coverage
#from code.functions import select_objects, requestData, trans4mEPSG, requestDataProcess, get_grid_indices_for_axis_X, get_grid_indices_for_axis_Y
from code.UserCred import saveCredentials
from code.GetLayers_Window import getLayers
from code.functions import *
import re
import os

# 1- Provide your user Credentials for fairicube.rasdaman.org if you have not saved them in an .env file yet
#saveCredentials("/path/to/your/envfile/.env")
#saveCredentials("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/.env")

# 2- Request some info about alyers available to select layers useful for your analysis, when providing a filepath, you can save this information as .csv
layer_info=getLayers(savepath="NONE")

# 3- Apply Coverage onto the *rasdaman* layers and make them ready to select
x=Coverage(layer_info)
x.getBoundary()

# 4- Filter out the samples that are covered by >1  of your selected layers 
#x.getSamples("path/to/samplefile/with/geoinformation/file.csv")

x.getSamples("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/dest_v2.samps_3May2024.csv")
x.getSamples("/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/samples_europe_pass.csv")
samplescovered=x.samples
#samplescorr={'AT_Kar_See_1_2016-08-01': ('46.8136889', '13.50794792'), 'AT_Nie_Mau_1_2015-07-20': ('48.375', '15.56'), 'AT_Nie_Mau_1_2015-10-19': ('48.375', '15.56'), 'AT_Wie_Gro_1_2012-08-03': ('48.2', '16.37'), 'AT_Wie_Gro_1_2012-10-20': ('48.2', '16.37')}


samplescorr=remove_partial_dates(samplescovered)
    #####

# 5- Create a list of layers by selecting them 
layers_to_analyze=select_objects("automatic",x)
info=layer_info[0]

layerlist=[]
for layer in layers_to_analyze:
    for i in range(0,len(x._data)):
        entries_list=list(x._data[i])
        try:
            index=entries_list.index(layer)
            layerlist.append(entries_list)
        except:
            continue 
    
out="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/wormpicked_data/rasdaman_data_test2.csv"
outstats="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/wormpicked_data/rasdaman_data_all_soilmoisture.csvStats.csv"
logpath="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/wormpicked_data/rasdaman_data_test2.log2.log"

# 6- Request the data for the SAMPLES for all chosen layers and save them as .csv
requestDataWGS(info,layerlist,samplescorr,out, logpath, offset=0, approximate=True)

#requestDataWGS(info,layerlist,samplescorr,"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/TestNALoop.csv", approximate="TRUE")

#reactivate the sys out
sys.stdout = sys.__stdout__
 
import pandas as pd
#data check 
summarize_csv(out,outstats)

##write Null values to list

data_list= layer_info_2811
#with open("NullValues_05_12_24.csv", "w") as file:
#    for item in data_list:
#        try:
#            # Get elements from index 12 and 13
#            entry_12 = item[12]  # Expected to be a list like ['a', 'b', 'c']
#            entry_13 = item[13]  # Expected to be a list like [1, 2, 3]
#            ln = item[0]
#            # Write each pair from entry_12 and entry_13 in a new row
#            for e12, e13 in zip(entry_12, entry_13):
#                file.write(f"{ln}_{e12},{e13}\n")
#                #print(f"{ln}_{e12},{e13}\n")
#        except IndexError:
#            # Handle the case where index 12 or 13 does not exist
#            print(f"Error: One of the entries at index 12 or 13 is missing in item: {item}")

with open("NullValues_10_12_24.csv", "w") as file:
    for item in data_list:
        try:
        # Get elements from index 12 and 13
            entry_12 = item[12]  # Expected to be a list like ['a', 'b', 'c']
            entry_13 = item[13]  # Expected to be a list like [1, 2, 3]
            ln = item[0]
            print(entry_12, entry_13)
        #print(entry_13)
        # Check if entry_13 is longer than entry_12 by 1 element, and if that element is 'NA'
            if len(entry_13) != len(entry_12):
                print(entry_13, entry_12) # Remove the 'NA' from entry_13
        except IndexError:
            # Handle the case where index 12 or 13 does not exist
            print(f"Error: One of the entries at index 12 or 13 is missing in item: {item}")