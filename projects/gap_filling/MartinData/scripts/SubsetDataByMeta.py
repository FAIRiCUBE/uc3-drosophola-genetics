import sys
import gzip
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--meta", dest="meta", help="Input file")
parser.add_option("--out", dest="out", help="Output file prefix")
parser.add_option("--criteria", dest="cri", help="Filtering criteria")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y


Criteria=options.cri.split("=")[0]
Arguments = options.cri.split("=")[1].split(",")

C=0
Keep=[]
O1=open(options.out+".meta","wt")
for l in load_data(options.meta):
    if l.startswith("sample"):
        header=l.rstrip().split()
        continue
    DATA = dict(zip(header, l.rstrip().split()))
    if DATA[Criteria] in Arguments:
        Keep.append(C)
        O1.write(l)
    C+=1

O2 = gzip.open(options.out+".af.gz", "wt")
for l in load_data(options.IN):
    a=l.rstrip().split()
    PL=a[:3]
    pops=a[3:]
    KL=[pops[x] for x in Keep]
    TEST=[sum([int(z) for z in x]) for x in zip(*[y.split(",") for y in KL])]
    if 0 in TEST:
        continue
    print("\t".join(PL)+"\t"+"\t".join(KL))



        
