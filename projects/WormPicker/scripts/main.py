import subprocess
import os
from scripts import module_crs_converter


os.chdir("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker")
subprocess.call(['python3', '/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/scripts/GetLayers.py'])

#output=$(python3 scripts/GetBoundary_New.py --LayerInfoFile output/layer_info_WCS.csv )

command = 'python3 scripts/GetBoundary_New.py --LayerInfoFile output/layer_info_WCS.csv'

# Run the command and capture the output
result = subprocess.check_output(command, shell=True)

# Decode the byte string result to a UTF-8 string
output = result.decode('utf-8')[:-1]

#command = f"python3 scripts/FilterSamples.py --source {source_file} --boundary {boundary}"
command = f"python3 /media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/scripts/FilterSamples.py --source /media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/data/dest_v2.samps_25Feb2023.csv --boundary {output}"

# Execute the command and capture the output
result = subprocess.run(command, shell=True, capture_output=True, text=True)

command= f"python3 /media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/scripts/GetDataNew.py --source/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/projects/WormPicker/output/coveredsamples_wcs.csv "
subprocess.run(command, shell=True)