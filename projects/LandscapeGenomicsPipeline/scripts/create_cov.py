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
parser.add_option("--input", dest="IN", help="Input file for variable reading")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--samples", dest="SAMP", help="samplename file where all the samplenames are listed")

(options, args) = parser.parse_args()
parser.add_option_group(group)

samples=[]
with open(options.IN, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        #print(line)
        samples.extend(line)


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
#indices_to_select= get_numeric_columns("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/metadata.csv")

#data_transposed = np.transpose(data)
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
#samps=filter_samples("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/metadata.csv", samples)
print(samps)
#data = [list(itemgetter(*indices_to_select)(row)) for i, row in enumerate(samps)]
data = [list(itemgetter(*indices_to_select)(row)) for i, row in enumerate(samps) if i > 0]

data_transposed = np.transpose(data)

with open(options.OUT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_transposed)

#info file
def get_colnames(x):
    with open(x, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        numeric_cols = []
        for i, row in enumerate(reader):
            for j, val in enumerate(row):
                if i == 0:
                    try:
                        float(val)
                        numeric_cols.append(j)
                    except ValueError:
                        pass
        numcols=[]
        for k in numeric_cols:
            numcols.append(header[k])
        return numcols

colnames = get_colnames(options.IN)

def get_colnames_write(x, bn):
    with open(x, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        numeric_cols = []
        for i, row in enumerate(reader):
            for j, val in enumerate(row):
                if i == 0:
                    try:
                        #float(val)
                        #print(header[j])
                        #print(data_transposed[j-1])
                        path= bn + "_" + header[j] + ".csv"
                        #print(path)
                        with open(path, 'w', newline='') as f:
                            array_as_string = ' '.join(map(str, data_transposed[j-1]))
                            f.write(array_as_string)
                    except ValueError:
                        pass
        return 0

with open(options.OUT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_transposed)


base_name, ext = os.path.splitext(options.OUT)
info_base_name = base_name + '.covariate.info'
options.OUT = info_base_name + ext
get_colnames_write(options.IN, info_base_name, )

with open(options.OUT, 'w') as f:
    print(*colnames, file=f)


