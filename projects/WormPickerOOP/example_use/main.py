from code.functions import Coverage
from code.module_crs_converter import trans4mEPSG
from code.Objects import Coverage
from code.functions import select_objects, requestData, trans4mEPSG
from code.UserCred import saveCredentials
from code.GetLayers import getLayers

# 1- Provide your user Credentials for fairicube.rasdaman.org
saveCredentials("/path/to/your/envfile/.env")
saveCredentials("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/.env4")

# 2- Request some info about alyers available to select layers useful for your analysis, when providing a filepath, you can save this information as .csv
layer_info=getLayers(savepath="NONE")

# 3- Apply Coverage onto the *rasdaman* layers and make them ready to select
x=Coverage(layer_info)
x.getBoundary()

# 4- Filter out the samples that are covered by >1  of your selected layers 
#x.getSamples("path/to/samplefile/with/geoinformation/file.csv")

x.getSamples("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv")
samplescovered=x.samples

# 5- Create a list of layers by selecting them 
layers_to_analyze=select_objects("manual",x)

layerlist=[]
for layer in layers_to_analyze:
    for i in range(0,len(x._data)):
        entries_list=list(x._data[i])
        try:
            index=entries_list.index(layer)
            layerlist.append(entries_list)
        except:
            continue

# 6- Request the data for the SAMPLES for all chosen layers and save them as .csv
requestData(layerlist,samplescovered, "/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/longtimetry.csv")




