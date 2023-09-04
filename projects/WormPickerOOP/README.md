
# Package Information

Usage: This package is mainly designed to assist download of rasterized data for scientific use. It was developed for the UC3 Project of FAIRiCUBE. The Wormpicker project is building the framework for  "Wormpicker" (UDF) of UC3.
Currently the OGC Service used is provided by rasdaman (https://ows.rasdaman.org/rasdaman/ows#/services), manual switch to another server host will be available but is not implemented yet. The FAIRiCUBE specific OGC is accessable via: https://fairicube.rasdaman.com/rasdaman/ows#/services 

## code 
All of the functions are located under the code.

## example use
How to use the package best is demonstrated in the example_use and a main.py.

# Asks for UserCredentials 
That are required for fairicube.rasdaman.org, will be written to an env file (path manually).

# Request Data
Please take into consideration, that at this time rasdaman layers do not provide a big variety of data in terms of time axis and therefore are not selected in time.
