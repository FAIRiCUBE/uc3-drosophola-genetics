import csv
import subprocess
import sys
from optparse import OptionParser, OptionGroup
import gzip
import os
import argparse
import ast


#########################################################   HELP   #########################################################################
usage="python %prog --source input.csv "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
""")
#########################################################   CODE   #########################################################################

parser.add_option("--LayerInfoFile", dest="INFO", help="The Layer Info file.")
parser.add_option("-l", "--list", dest="my_list", action="append", help="Usage: python my_script.py -l Object1 -l Object2 -l Object3")

parser.add_option_group(group)
(options, args) = parser.parse_args()

layer_names= options.my_list

#print(layer_names)

# Open the CSV file and read its contents
with open(options.INFO, "r") as f:
    reader = csv.reader(f)
    next(reader)  # Skip the header row
    rows = list(reader)

# Extract the MinX, MaxX, MinY, and MaxY values for the specified layers
min_x_values = []
max_x_values = []
min_y_values = []
max_y_values = []

for layer_name in layer_names:
    layer_rows = [row for row in rows if row[0] == layer_name]
    if len(layer_rows) > 0:
        layer_min_x_values = [float(row[1]) for row in layer_rows]
        layer_max_x_values = [float(row[2]) for row in layer_rows]
        layer_min_y_values = [float(row[3]) for row in layer_rows]
        layer_max_y_values = [float(row[4]) for row in layer_rows]
        min_x_values.append(min(layer_min_x_values))
        max_x_values.append(max(layer_max_x_values))
        min_y_values.append(min(layer_min_y_values))
        max_y_values.append(max(layer_max_y_values))

# Calculate the common boundary box
common_min_x = min(min_x_values)
common_max_x = max(max_x_values)
common_min_y = min(min_y_values)
common_max_y = max(max_y_values)

# Print the common boundary box
print(f"{common_min_x}:{common_max_x}:{common_min_y}:{common_max_y}")
