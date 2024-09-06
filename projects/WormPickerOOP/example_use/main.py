from code.module_crs_converter import trans4mEPSG
from code.Objects import Coverage
#from code.functions import select_objects, requestData, trans4mEPSG, requestDataProcess, get_grid_indices_for_axis_X, get_grid_indices_for_axis_Y
from code.UserCred import saveCredentials
from code.GetLayers_Window import getLayers
from code.functions import *
import re

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
x.getSamples("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/EuropeSamples.csv")
samplescovered=x.samples


samplescorr=remove_partial_dates(samplescovered)
    #####

# 5- Create a list of layers by selecting them 
layers_to_analyze=select_objects("manual",x)
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
    
out="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/rasdaman_data.csv"
outstats="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/RasdamanDataAllSept24._Stats.csv"
logpath="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/rasdaman_data.log"

# 6- Request the data for the SAMPLES for all chosen layers and save them as .csv
requestDataWGS(info,layerlist,samplescorr,out, logpath, offset=0, approximate=True)

#requestDataWGS(info,layerlist,samplescorr,"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/TestNALoop.csv", approximate="TRUE")


#reactivate the sys out
sys.stdout = sys.__stdout__
 
import pandas as pd
#data check 
summarize_csv(out,outstats)


