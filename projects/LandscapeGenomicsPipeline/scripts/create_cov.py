import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import csv
from csv import reader, writer
import numpy as np
from operator import itemgetter
import os


#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################
parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--samples", dest="SAMP", help="samplename file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

samples=[]
with open(options.SAMP, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        #print(line)
        samples.extend(line)


##with open(options.IN, 'r', encoding='latin-1') as f:
##    reader = csv.reader(f, delimiter=',')
##    data = [row[29:] for i, row in enumerate(reader) if i > 0 and row[0] in samples]
##
##data_transposed = np.transpose(data)
##
##with open(options.OUT, 'w', newline='') as f:
##    writer = csv.writer(f, delimiter=' ')
##    writer.writerows(data_transposed)
##
###new version going from metafile


def get_numeric_columns(x):
    with open(x, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        #print(header)
        numeric_cols = []
        for i, row in enumerate(reader):
            #print(i)
            for j, val in enumerate(row):
                #print(j)
                #print(val)
                if i == 0:
                    #print(row)
                    #print(val)
                    try:
                        float(val)
                        numeric_cols.append(j)
                    except ValueError:
                        pass
                #elif j in numeric_cols:
                    #try:
                    #    float(val)
                    #except ValueError:
                    #    numeric_cols.remove(j)
        numcols=[]
        for k in numeric_cols:
            numcols.append(header[k])
        return numeric_cols

indices_to_select= get_numeric_columns(options.IN)


def filter_samples(meta, samples):
    filtered_rows= []
    with open (meta, 'r') as file:
        reader= csv.reader(file)
        for row in reader:
            if row and row[0] in samples:
                filtered_rows.append(row)
    return filtered_rows

print(samples)
samps = filter_samples(options.IN, samples)
print(samps)
data = [list(itemgetter(*indices_to_select)(row)) for i, row in enumerate(samps)]

#with open(options.IN, 'r', encoding='latin-1') as f:
#    reader = csv.reader(f, delimiter=',')
#    data = [list(itemgetter(*indices_to_select)(row)) for i, row in enumerate(reader) if i > 0]

data_transposed = np.transpose(data)

with open(options.OUT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_transposed)

#info file
def get_colnames(x):
    with open(x, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        #print(header)
        numeric_cols = []
        for i, row in enumerate(reader):
            #print(i)
            for j, val in enumerate(row):
                #print(j)
                #print(val)
                if i == 0:
                    #print(row)
                    #print(val)
                    try:
                        float(val)
                        numeric_cols.append(j)
                    except ValueError:
                        pass
                #elif j in numeric_cols:
                    #try:
                    #    float(val)
                    #except ValueError:
                    #    numeric_cols.remove(j)
        numcols=[]
        for k in numeric_cols:
            numcols.append(header[k])
        return numcols

colnames = get_colnames(options.IN)

base_name, ext = os.path.splitext(options.OUT)
info_base_name = base_name + '.used.info'
options.OUT = info_base_name + ext

with open(options.OUT, 'w') as f:
    print(*colnames, file=f)


