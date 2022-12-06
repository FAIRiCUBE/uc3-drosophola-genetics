import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip

    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


for l in load_data(options.IN):
    a = l.rstrip().split()
    if l.startswith("##"):
        continue
    if l.startswith("#"):
        header = a[9:]
        print("Chr\tPos\t" + "\t".join(header))
        continue
    pops = a[9:]
    format = a[8].split(":")
    AFs = []
    for i in pops:
        if "./." in i:
            AFs.append("NA")
            continue
        P = dict(zip(format, i.split(":")))
        AFs.append(str(round(float(P["AD"]) / float(P["DP"]), 3)))
    if sum([float(x) for x in AFs if x != "NA"]) == 0:
        continue
    print(a[0] + "\t" + a[1] + "\t" + "\t".join(AFs))
