from code.functions import getSamples
from code.module_crs_converter import trans4mEPSG
from code.Objects import Coverage
from code.functions import show_selected_options, select_objects, requestData

from code.UserCred import saveCredentials
from code.GetLayers import getLayers

# 1- Provide your user Credentials for fairicube.rasdaman.org
saveCredentials("/path/to/your/envfile/.env")
layer_info=getLayers(savepath="NONE")
x=Coverage(layer_info)
x.getBoundary()
x.getSamples("path/to/samplefile/with/geoinformation/file.csv")
samplescovered=x.samples

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

requestData(layerlist,samplescovered, "/result/filepath")





