# Purpose
Download environmental data from the Rasdaman GmbH database for predefined point coordinates for a specific date if available.



# Web Application

If you want to use the WP Tool via WebApp, please visit: 
https://wormpicker.nilu.no/
Required input file is a .csv file with columns "sampleId, lat, long, date".



# Package Information

Usage: This package is mainly designed to assist download of rasterized data for scientific use. It was developed for the UC3 Project of FAIRiCUBE.
Currently the OGC Service used is provided by rasdaman (https://ows.rasdaman.org/rasdaman/ows#/services), manual switch to another server host will be available but is not implemented yet. The FAIRiCUBE specific OGC is accessable via: https://fairicube.rasdaman.com/rasdaman/ows#/services 

## Query Cube
While many applications are leveraging the gridded datasets becoming available from diverse sources such as Copernicus services, the European Environment Agency or Eurostat, in most cases, these focus on the spatiotemporal extent, with analyses being performed across multiple grid cells. In addition, the wealth of data provided by these sources can also be valuable for applications focusing on point-based locations. However, at present, it is difficult and resource intensive to extract values for the relevant point-based locations from the gridded resources.

The FAIRiCUBE "Query Cube" aims to fill this gap, enabling users with little experience in working with geospatial data to access relevant values from gridded resources for the point locations they’re interested in. All the user must do is to provide a coordinate pair, and indicate for which of the available gridded resources they’d require values. The Worm Picker accesses the complete FAIRiCUBE at the spatial location specified, and extracts a vector containing the requested information.

<img src="https://github.com/FAIRiCUBE/uc3-drosophola-genetics/assets/11915304/a5357f3b-2ebf-4cb4-8f41-83ec4c67644b" width="281" height="430" />

Query Cube is being developed utilizing the OGC Web Coverage Service (WCS), enabling dynamic access to the data being stored on the FAIRiCUBE rasdaman instance. This software will be made available in different forms:
- Web GUI: enables the user to interactively specify both the spatial location and the relevant gridded resources for which to return data;
- Python Library: enables integration of the Worm Picker functionality within python scripts;

This code has been developed within UC3 of FAIRiCUBE the Project. 

Currently the OGC WCS Service used is provided by rasdaman (https://ows.rasdaman.org/rasdaman/ows#/services).
The FAIRiCUBE specific WCS is accessable via: https://fairicube.rasdaman.com/rasdaman/ows#/services 

The main workflow established is executed by several python scripts manually, documented compact in a [main file](example_use/main.py) with execution code, further explanation and backround can be read under [Workflow](#workflow). 



## Package Structure

**Requirements**

- *WormpickerEnv.yml* - To set up an environment with all the required installations. 

**Directories**
- *code*: Stores all the python files containing functions and objects that are required for the Wormpicker to work. 
- *example_use*: Contains a main.py file, that calls all the functions required to download data for geo referenced samples as well as example data file called *samplesfile.csv*.

# Functions

All the necessary functions can be found in the [code](projects/QueryCube/code) directory.


# Example Use


## Running Query Cube 

The program can be run from the terminal with the following command line:


```
 python3 main.py \
    -outdir example_use/example_results \
    -samples samplesfile.csv \
    -username "rasdaman_username" \
    -password "rasdaman_password" 
 
```

The following parameters can be used with the script. 

- Rasdaman credentials: 
    - Username
    - Password
    - Service Endpoint

    If not provided manually , these parameters can be read from an env file as well.

- samples: A csv file carrying information on sampleId, lat, long and date. This file is mandatory. 
- outdir: If the result should be stored as csv provide the path to the output directory. 


Additional functions in the modules can be used to run custom workflows of Query Cube.


# Output
The desired output when using this package is a .csv file, that is comprised of the grid cell data for each desired layer (columns) for all the geo referenced samples (rows).
Please provide the path to you desired output folder in the "RequestData" step. 

The purpose of this structure origiated by the ecessity to use the .csv file in a downstream pipeline for Environmental Association Analysis, and further to facilitate intersecting the grid cell values with other metadata to perform regression analysis. 



