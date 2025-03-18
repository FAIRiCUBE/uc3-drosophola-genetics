import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--GODescriptions", dest="GD", help="Input file")
parser.add_option("--GOAnnotations", dest="GA", help="Input file")

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


GOs = d(str)

for l in load_data(options.GD):
    if l.startswith("id: "):
        ID = l.rstrip().split("id: ")[1]
    if l.startswith("name: "):
        NAME = l.rstrip().split("name: ")[1]
        GOs[ID] = NAME
    if l.startswith("alt_id: "):
        ID = l.rstrip().split("alt_id: ")[1]
        GOs[ID] = NAME

GENEs = d(list)
for l in load_data(options.GA):
    if l.startswith("!"):
        continue
    a = l.rstrip().split("\t")
    GENEs[a[4]].append(a[1])

for k, v in GOs.items():
    if k not in GENEs:
        continue
    print(k, v, " ".join(GENEs[k]), sep="\t")
