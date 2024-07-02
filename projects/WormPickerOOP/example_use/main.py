from code.functions import Coverage
from code.module_crs_converter import trans4mEPSG
from code.Objects import Coverage
from code.functions import select_objects, requestData, trans4mEPSG, requestDataProcess, get_grid_indices_for_axis_X, get_grid_indices_for_axis_Y
from code.UserCred import saveCredentials
from code.GetLayers import getLayers
import code.functions

# 1- Provide your user Credentials for fairicube.rasdaman.org
saveCredentials("/path/to/your/envfile/.env")
saveCredentials("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/.env")

# 2- Request some info about alyers available to select layers useful for your analysis, when providing a filepath, you can save this information as .csv
layer_info=getLayers(savepath="NONE")

# 3- Apply Coverage onto the *rasdaman* layers and make them ready to select
x=Coverage(layer_info)
x.getBoundary()

# 4- Filter out the samples that are covered by >1  of your selected layers 
#x.getSamples("path/to/samplefile/with/geoinformation/file.csv")

x.getSamples("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/dest_v2.samps_3May2024.csv")
samplescovered=x.samples

import re
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
    
out="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/example_use/WGSBasedResult.csv"

# 6- Request the data for the SAMPLES for all chosen layers and save them as .csv
requestDataWGS(info,layerlist,samplescorr,out)




