
# Project Eggdance

## Summary 
For the WormPicker Deliverable of FAIRiCUBE, point specific layer data for geospatial points of interest is used, when operating on the rasdaman OGC service endpoint (https://fairicube.rasdaman.com/rasdaman/ows#/services). There are two options of retrieving data via request, the first one being „slicing“ and the second one „trimming“.

Slicing is for a sample data retrievel sufficient, assuming that each point adressed is retrieving a data value. In some cases and also furhter development of data retrieval, trimming would be essential, to establish an area of a desired size around the geospatial point desired (e.g. a 100m radius around this point) to get a realisitic estimate of the grid values around it. 
The basic assumtion was, that taking a grid point and trimming it will work with an fold-increase oft he grid resolution per se. Starting with information about (depicted via browser interface here) a layer, we can assess that the CRS of this particular layer is EPSG/0/3035, the layer covers an area from Y (Latitude in WGS 84) 900000m to 5500000m and on X (Longitude in WGS 84) from 900000m to 7400000m and the grid cell resoltion is (100 * -100 metre), see the layer_info.csv file [here](/projects/WormPicker/output/new_layer_info_WCS.csv). 

Assumption: A grid cell covers an area of 100m * 100m.


Slicing at the exact sample location point (in this example (4588860.382936129 at X; 2639692.4010280264 at Y) does deliver one value, the same is true for other samples with other geospatial coordinates and other layers. For  this specific one the result is the Null-Value  {-128}. The question arises, since tried for many other points, if the whole layer is actually covered in NullValue containing grid cells. The Slicing operation in the browser interface is depicted below. For using the python „request“ package, a substring oft he request would include „&subset=Y({},{})&subset=X({},{})" with corresponding values in the curly brackets. Slicing would use the exact same value twice, whereas trimming would need a higher and lower value for each X and Y.

**Does trimming deliver a useful result?**
This question hast o be subcategorized into several single questions and approaches to understand the mechanistics of layer/grid data retrieval.

*Approaches to solve this problem:*
Thie third representation in the row above is trying to depict the problem of what in the projects context would be „directionality“, meaning which direction the grid gets sliced into. Going „positive“ meters on both X and Y woudl represent a north-eastern slicing when using the sample point as reference. Since the projection dynamics of the CRS are not fully clarified yet, the directionality will now be referred to as quadrant, since the four combination on positive and negative on two dimensions like X and Y resemble the mathematical concept of quadrants (Q1, Q2, Q3, Q4 = (1,1), (-1,1), (-1,-1), (1,-1)).

**How much trimming is useful/possible and into what is a useful directionality?**
„Question 1“ was tested by plotting how much cells are retrieved by the rasdaman query when different „directions“ are used to trim the grid (area of each square = number of values). The assumption was already explained, that a 100 m x 100 m trimming represents one gridcell (resolution 100m), therefore 200m x 200m should deliver 4 grid cells (in this particular layer). This was tested for various samples within the layer.  1x and 2x increase of the trimming deliver equal results for all combinations (see below). Fort he 2x resolution, it is obvoious that the query delivers four values, this pattern is consistent for all natural steps of resolution, meaning that integer numbers are delivering equally distrubuted results. This was leading tot he next question.

**What if the fold increase of the „radius“ is not an integer, but a decimal value/float?**
When considering a radius in terms of biological/ecological questions, can a defined size be used (e. g. 150 m around the sampling point)? This would represent 1,5 times the resolution oft he grid, how ist he distribution of values in that case? 
The preliminary answer to this question is essential. Depending on where in the „original“ grid cell the exact point is located, going a certain radius above 1x resolution will deliver different results in respect to the direction, where the trimming is happening to. This means the following:
*Trimming the layer from sampling point (AT_Kar_see) +150m on the X axis and +150m in respect to the sampling point on the Y axis will deliver one value, whereas trimming on -150m on the X axis and -150m in respect to the point will deliver four values. This was tried for different samples as well and delivered different results in terms of „directionality“, meaning which direction was delivering most values. In the Austrian sample point we see „the oppsosite behaviour“ of the querying in comparison to a sample from spain or a serbian geopoint.* 

This behaviour was tested for different fine-gridded versions of resultion, and preliminary assumptions are, that after x.75 resolution the values are „rounded upwards“. One approach to determine where this differences in directionality are coming from was, wether or not the position in the grid cell, defined by the modulo value of X and Y (X%100 and Y%100) is causing these differences.

To determine, representative „Testpoints“ for certain combinations of grid-cell locations were established and tested (defined in quarter-steps, T1-T9) and this set was repeated for four different locations in the layer (lower left value, lower right value, upper left value, upper right value). They were changed manually, the script contains the coding for one of these four sets (https://github.com/FAIRiCUBE/uc3-drosophola-genetics/blob/NHM/projects/EggTest/EggTest.py#L64-L72).
This was done to get a better insight, if the value bias is changing depending on the grid position or if it stays the same across the whole layer. 
The lower points are not sitting on the border oft he layer, in fact the are 500m away from the border but the graph is plotting the points directly on the blayer border. Each T1-T9 corridnates are described below. 
