import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import operator

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--window", dest="WI", help="Input file")
parser.add_option("--chromatin", dest="CHROM", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


ChromatinD = d(lambda: d(list))
ChromatinL = d(lambda: d(list))
Window = d(lambda: d(list))
Chromosome = d(list)
WI = int(options.WI)
ops = {">": operator.gt, "<": operator.lt}  # etc.


for l in load_data(options.CHROM):
    if l.startswith("Chrom"):
        continue
    a = l.rstrip().split()
    ChromatinD[a[0]][a[2]].append(a[1])

for l in load_data(options.IN):
    if l.startswith("#"):
        continue
    a = l.rstrip().split()
    TEST = 0
    for Type, v1 in ChromatinD[a[0]].items():
        for Pos in v1:
            if ops[Pos[0]](int(a[1]), float(Pos[1:])*1000000):
                # print(Pos[0], int(a[1]), float(Pos[1:])*1000000)
                ChromatinL[Type][a[0]].append(int(a[1]))
                TEST = 1
    if TEST == 0:
        ChromatinL["Euchromatin"][a[0]].append(int(a[1]))
    Chromosome[a[0]].append(int(a[1]))
    Window[a[0]+"_"+str(int(round(float(a[1])/WI, 0)*WI))][a[0]].append(a[1])

OUT1 = open(options.OUT+"_windows.csv", "w")
OUT1.write("Type,Chrom,Pos\n")
for Type, v in sorted(Window.items()):
    for Chrom, v1 in sorted(v.items()):
        for Pos in sorted(v1):
            OUT1.write(Type+","+Chrom+","+str(Pos)+"\n")
OUT1.close()
OUT2 = open(options.OUT+"_chromatin.csv", "w")
OUT2.write("Type,Chrom,Pos\n")
for Type, v in sorted(ChromatinL.items()):
    for Chrom, v1 in sorted(v.items()):
        for Pos in sorted(v1):
            OUT2.write(Type+","+Chrom+","+str(Pos)+"\n")
OUT2.close()
OUT3 = open(options.OUT+"_chromosomes.csv", "w")
OUT3.write("Type,Chrom,Pos\n")
for Chrom, v1 in sorted(Chromosome.items()):
    for Pos in sorted(v1):
        OUT3.write(Chrom+","+Chrom+","+str(Pos)+"\n")
OUT3.close()
