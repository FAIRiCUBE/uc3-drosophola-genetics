B:
In the attachment is the python script file which represents the formula (with EPSG:4326 coordinates) how petascope translates input geo-subsets from WCS GetCoverage request to grid indices (zero-based index) to query from imported rasdaman collection. Note that there is rounding from decimal numbers to integer grid indices.
You can adjust it with the geo extents and geo resolutions of X Y axes from coverage which you want to test to get the grid indices.

import os
import math

# nrow = 10, col = 10, ncell = 100
# resx = 36, resy = 18
# extent = -180, 180, -90, 90 (xmin, xmax, ymin, ymax)

# It works like in R trunc(float number)
def shiftToInteger(float_number):
    rounded_number = math.ceil(float_number)
    if float_number < 0:
        rounded_number = math.floor(float_number)
        
    tmp = abs(float_number) + 0.000000000000000222
    
    if tmp >= abs(rounded_number):
        return rounded_number

    return int(float_number)        

# geo extent in EPSG:4326
xmin = -180
xmax = 180
ymin = -90
ymax = 90

# resolutions of Long (X) and Lat (Y) axes
xres = 36
yres = 18

"""
def colFromX(x):
    colnr = int( (x - xmin) / xres ) + 1
    return colnr

def rowFromY(y):
    rownr = 1 + ( int( (ymax - y) / yres ) )
    return rownr
"""

def colFromX(x):
    colnr = shiftToInteger( (x - xmin) / xres ) + 1
    return colnr

def rowFromY(y):
    rownr = 1 + ( shiftToInteger( (ymax - y) / yres ) )
    return rownr
    

def cal(subset_xmin, subset_ymin, subset_xmax, subset_ymax):
    col1 = colFromX(subset_xmin + 0.5 * xres) 
    col2 = colFromX(subset_xmax - 0.5 * xres)
    row1 = rowFromY(subset_ymax - 0.5 * yres) 
    row2 = rowFromY(subset_ymin + 0.5 * yres)
    
    print("grid indices for X:", col1 - 1, col2 - 1, "\ngrid indices for Y:", row1 - 1, row2 - 1)


# subset by geo extents: minX, minY, maxX, maxY        
cal(-180, -90, 180, 90)



S: My understanding would be the following: (grid indices x: 0 1; grid indices y: 0 1  1 grid cell queried?) and e.g. (grid indices x: 1 2, grid indices y: 1 2  1 grid cell queried  -  whereas grid indices x: 0 2, grid indices y: 0 2  4 grid cells queried?) 
And if so, would you advise to use the resolution of a grid as estimator to use “trimming” to retrieve point data or another value (like 0.5*resolution) to keep the within-grid-indices-distance as similar as possible for multiple requested coordinates? 

B: A request with fixed geo subsettings on X and Y axes (e.g. via WCS GetCoverage request) cannot return multiple results.
In the python script which I sent to you, it has the formula to return the corresponding grid indices (these values are deterministic from the formula).
In rasdaman, grid index starts with 0 and a trimming with grid indices:
- x[0:0], y[0:0] -> 1 grid cell
- x[0:1], y[0:0] -> 2 grid cells
- x[0:0], y[0:1] -> 2 grid cells
- x[0:1], y[0:1] -> 4 grid cells
...
- x[0:2], y[0:2] -> 9 grid cells ([0:2] means [0,1,2])


S: We wanted to query the given subsets, why does a subset of size 5x5 or 15x15 provide two values/indices, while a subset of 10x10 only provides one?
Maybe I oversaw and you already wrote, but just to be sure what’s going on: Where does the 0.000000000000000222 in the tmp=… formula come from?

Output of cal()-functions:
>>> cal(15,15,15,15)
grid indices for X: 20 19 
grid indices for Y: 8 7
>>> cal(10,10,15,15)
grid indices for X: 19 19 
grid indices for Y: 8 7
>>> cal(10,10,20,20)
grid indices for X: 19 19 
grid indices for Y: 7 7
>>> cal(10,10,25,25)
grid indices for X: 19 20 
grid indices for Y: 7 7

def shiftToInteger_print(float_number):
     print("Your input is:", float_number)
     rounded_number = math.ceil(float_number)
     print("I have rounded to:", rounded_number)
     if float_number < 0:
         rounded_number = math.floor(float_number)
         print("Input Number is less than 0, now rounded number is", rounded_number)
     tmp = abs(float_number) + 0.000000000000000222
     print("Tmp is:", tmp)
     if tmp >= abs(rounded_number):
         print("Tmp >= abs(rn)", rounded_number)
         return rounded_number
     print(int(float_number))
 
>>> shiftToInteger_print(0.5)
Your input is: 0.5
I have rounded to: 1
Tmp is: 0.5000000000000002
0



B:
Transform from spatial geo extents to grid extents are not easy because people has different perspectives about the output grid domains (+ or - 1 grid pixel).
This constant value 0.000000000000000222 comes from my tests in R language for rounding float number to integer. I tried to resemble R trunc() function:
> trunc(4.0 - 0.000000000000000222)
[1] 4
> trunc(4.0 - 0.00000000000000222)
[1] 3
I can tell you only from my experience with tries and errors, there are two cases for these geo-to-grid conversion which I updated in the script at the bottom (here grid lower bound starts with 0 to be simplified):
- First case, (subset geo max - subset geo min) > (geo res / 2), then use formula for bigger than half grid pixel
- Second case, (subset geo max - subset geo min) < (geo res / 2), then use formula for smaller than half grid pixel. This case is special, because gdal_translate for example does not support to subset less than half grid pixel.
For example: test.tiff has 10 x 10 grid pixels and geo extents with Long(-180:180) and Lat(-90:90)
gdal_translate -projwin -180 90 -170 85 test.tiff subset.tiff
Input file size is 10, 10
ERROR 1: Error: Computed -srcwin 0 0 0 0 has negative width and/or height.

The updated script with your examples and xres=10 degrees for Lon axis and yres=-10 degrees for Lat axis is here:
import os
import math

# nrow = 10, col = 10, ncell = 100
# resx = 36, resy = 18
# extent = -180, 180, -90, 90 (xmin, xmax, ymin, ymax)

# It works like in R trunc(float number)
def shiftToInteger(float_number):
    rounded_number = math.ceil(float_number)
    if float_number < 0:
        rounded_number = math.floor(float_number)
        
    tmp = abs(float_number) + 0.000000000000000222
    
    if tmp >= abs(rounded_number):
        return rounded_number

    return int(float_number)        

# geo extent in EPSG:4326
xmin = -180
xmax = 180
ymin = -90
ymax = 90

# resolutions of Long (X) and Lat (Y) axes
xres = 10

# NOTE: Y axis (negative resolution), grid origin starts from max geo bound
# (e.g. Lat axis, then pixel 0 is at 90 degree and pixel 9th at -90 degree with yres = -18)
yres = -10

"""
def colFromX(x):
    colnr = int( (x - xmin) / xres ) + 1
    return colnr

def rowFromY(y):
    rownr = 1 + ( int( (ymax - y) / yres ) )
    return rownr
"""

def colFromX(x):
    colnr = shiftToInteger( (x - xmin) / xres )
    return colnr

def rowFromY(y):
    rownr = ( shiftToInteger( (ymax - y) / abs(yres) ) )
    return rownr


def cal_less_than_half_pixel_X(subset_xmin, subset_xmax):
    grid_lower_bound_X = (subset_xmin - xmin) / xres
    grid_lower_bound_X = math.floor(grid_lower_bound_X)

    grid_upper_bound_X = (subset_xmax - xmin) / xres
    grid_upper_bound_X = math.ceil(grid_upper_bound_X) - 1

    return grid_lower_bound_X, grid_upper_bound_X

def cal_less_than_half_pixel_Y(subset_ymin, subset_ymax):
    grid_lower_bound_Y = (subset_ymax - ymax) / yres
    grid_lower_bound_Y = math.floor(grid_lower_bound_Y)

    grid_upper_bound_Y = (subset_ymin - ymax) / yres
    grid_upper_bound_Y = math.ceil(grid_upper_bound_Y) - 1

    return grid_lower_bound_Y, grid_upper_bound_Y


def cal_more_than_half_pixel_X(subset_xmin, subset_xmax):
    col1 = colFromX(subset_xmin + 0.5 * xres)
    col2 = colFromX(subset_xmax - 0.5 * xres)

    grid_lower_bound_X = col1
    grid_upper_bound_X = col2

    if grid_lower_bound_X > grid_upper_bound_X:
        grid_lower_bound_X = grid_upper_bound_X

    return grid_lower_bound_X, grid_upper_bound_X


def cal_more_than_half_pixel_Y(subset_ymin, subset_ymax):
    row1 = rowFromY(subset_ymax - 0.5 * abs(yres))
    row2 = rowFromY(subset_ymin + 0.5 * abs(yres))

    grid_lower_bound_Y = row1
    grid_upper_bound_Y = row2

    if grid_lower_bound_Y > grid_upper_bound_Y:
        grid_lower_bound_Y = grid_upper_bound_Y

    return grid_lower_bound_Y, grid_upper_bound_Y


def cal(subset_xmin, subset_xmax, subset_ymin, subset_ymax):

    if (subset_xmax - subset_xmin) < (xres / 2):
        grid_lower_bound_X, grid_upper_bound_X = cal_less_than_half_pixel_X(subset_xmin, subset_xmax)
    else:
        grid_lower_bound_X, grid_upper_bound_X = cal_more_than_half_pixel_X(subset_xmin, subset_xmax)


    if (subset_ymax - subset_ymin) < (abs(yres) / 2):
        grid_lower_bound_Y, grid_upper_bound_Y = cal_less_than_half_pixel_Y(subset_ymin, subset_ymax)
    else:
        grid_lower_bound_Y, grid_upper_bound_Y = cal_more_than_half_pixel_Y(subset_ymin, subset_ymax)

    print("Input geo bounds X: " + str(subset_xmin) + ":" + str(subset_xmax))
    print("Input geo bounds Y: " + str(subset_ymin) + ":" + str(subset_ymax))
    print("Calculated grid bounds for X: " + str(grid_lower_bound_X) + ":" + str(grid_upper_bound_X))
    print("Calculated grid bounds for Y: " + str(grid_lower_bound_Y) + ":" + str(grid_upper_bound_Y))


# subset by geo extents: minX, maxX, minY, maxY

cal(-180, 180, -90, 90)
print("----------------")

cal(15, 15, 15, 15)
print("----------------")

cal(10, 15, 10, 15)
print("----------------")

cal(10, 20, 10, 20)
print("----------------")

cal(10, 25, 10, 25)
print("----------------")

cal(155, 180, 46, 90)
print("----------------")


