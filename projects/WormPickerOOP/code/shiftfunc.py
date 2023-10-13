import os
import math

####???

### converting geographic coordinates to grid cell indices within a specified grid

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
yres = 10

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
cal(-180, -90, 180, 90)#
cal(-180,80,-170,90)
cal(0,0,0,0)
cal(10,10,20,20)
cal(10,10,25,25)


shiftToInteger(0.5)

# It works like in R trunc(float number)
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
    #return int(float_number)        