from argparse import ArgumentParser
from code.module_crs_converter import trans4mEPSG
from code.Objects import *
from code.UserCred import saveCredentials
from code.GetLayers_Window import getLayers
from code.functions import *
import re
import os

log_object = LogObject()
logger = log_object.logger

argparse = ArgumentParser()
argparse.add_argument("-outdir", help="Path to the outputfolder.")
argparse.add_argument("-samples", help='Samplesfile that carries information on ID, lat and long', required=True)
argparse.add_argument("-username", help='Username for Rasdaman')
argparse.add_argument("-password", help='Password for Rasdaman')
argparse.add_argument("-endpoint", help='Service Endpoint e.g. Rasdaman')
#argparse.add_argument("-samples", help='Samplesfile that carries information on ID, lat and long', required=True)

args = argparse.parse_args()
outdir = args.outdir
samples = args.samples
username = args.username
password = args.password
endpoint = args.endpoint


# A) Provide your user Credentials for fairicube.rasdaman.org if you have not saved them in an .env file yet
# B) Request some info about alyers available to select layers useful for your analysis, when providing a filepath, you can save this information as .csv
#layer_info=getLayers(savepath="NONE", loggerobject=logger)

layer_info = getLayers(savepath="NONE",
                       rasdaman_endpoint=endpoint,
                       rasdaman_password=password,
                       rasdaman_username=username,
                       loggerobject=logger)

# 3- Apply Coverage onto the *rasdaman* layers and make them ready to select
x=Coverage(layer_info)
x.getBoundary()
#x.getSamples("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/dest_v2.samps_3May2024.csv")
#x.getSamples("/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/samplesfile.csv")
x.getSamples(samples)
samplescovered=x.samples
#samplescorr={'AT_Kar_See_1_2016-08-01': ('46.8136889', '13.50794792'), 'AT_Nie_Mau_1_2015-07-20': ('48.375', '15.56'), 'AT_Nie_Mau_1_2015-10-19': ('48.375', '15.56'), 'AT_Wie_Gro_1_2012-08-03': ('48.2', '16.37'), 'AT_Wie_Gro_1_2012-10-20': ('48.2', '16.37')}
#samplescorr=remove_partial_dates(samplescovered)

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

if outdir:    
    outvar = outdir+"/WormpickerResult.csv"
    logpath = outdir+"/WormpickerResult.log"
else:
    outvar="NONE"

# Request the data for the SAMPLES for all chosen layers and save them as .csv
# Credentials taken from .env file, you can save result as object
# requestDataWGS(info,layerlist,samplescovered,"NONE", offset=0, approximate=True, loggerobject=logger)

# Credentials parsed from input, results can be saved as object with e.g.: result = requestDataWGS()....
requestDataWGS(info,layerlist,samplescovered,filepath=outvar, offset=0, approximate=True, rasdaman_password=password, rasdaman_username=username, rasdaman_endpoint=endpoint, loggerobject=logger)


with open(logpath, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    for log in log_object.get_logs():
        #print(log_object.get_logs()[log])
        row=[str(log)]
        writer.writerow(row)

#requestDataWGS(info,layerlist,samplescorr,"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/TestNALoop.csv", approximate="TRUE")

##reactivate the sys out
#sys.stdout = sys.__stdout__
# 
##data check 
#summarize_csv(out,outstats)
#
###write Null values to list
#
#data_list= layer_info
#
#with open("NullValues_10_12_24.csv", "w") as file:
#    for item in data_list:
#        try:
#        # Get elements from index 12 and 13
#            entry_12 = item[12]  # Expected to be a list like ['a', 'b', 'c']
#            entry_13 = item[13]  # Expected to be a list like [1, 2, 3]
#            ln = item[0]
#            print(entry_12, entry_13)
#        #print(entry_13)
#        # Check if entry_13 is longer than entry_12 by 1 element, and if that element is 'NA'
#            if len(entry_13) != len(entry_12):
#                print(entry_13, entry_12) # Remove the 'NA' from entry_13
#        except IndexError:
#            # Handle the case where index 12 or 13 does not exist
#            print(f"Error: One of the entries at index 12 or 13 is missing in item: {item}")


#follow()