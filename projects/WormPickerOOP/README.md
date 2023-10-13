
# Package Information

Usage: This package is mainly designed to assist download of rasterized data for scientific use. It was developed for the UC3 Project of FAIRiCUBE. The Wormpicker project is building the framework for  "Wormpicker" (UDF) of UC3.
Currently the OGC Service used is provided by rasdaman (https://ows.rasdaman.org/rasdaman/ows#/services), manual switch to another server host will be available but is not implemented yet. The FAIRiCUBE specific OGC is accessable via: https://fairicube.rasdaman.com/rasdaman/ows#/services 


## Package Structure

**Requirements**
*WormpickerEnv.yml* - To set up an environment with all the required installations. 

**Directories**
- *code*: Stores all the python files containing functions and objects that are required for the Wormpicker to work. 
- *example_use*: Contains a main.py file, that calls all the functions required to download data for geo referenced samples as well as example data file called *dest_v2.samps_25Feb2023.csv* (https://dest.bio/).


## Example use
How to use the package best is demonstrated in the example_use and a main.py.

1. 
2. 

# Functions

## Asks for UserCredentials 
That are required for fairicube.rasdaman.org, will be written to an env file (path manually).

## Request Data
Please take into consideration, that at this time rasdaman layers do not provide a big variety of data in terms of time axis and therefore are not selected in time.


# Unterstanding The Output
The desired output when using this package is a .csv file, that is comprised of the grid cell data for each desired layer (columns) for all the geo referenced samples (rows).
Please provide the path to you desired output folder in the "RequestData" step. 

The purpose of this structure origiated by the ecessity to use the .csv file in a downstream pipeline for Environmental Association Analysis, and further to facilitate intersecting the grid cell values with other metadata to perform regression analysis. 



