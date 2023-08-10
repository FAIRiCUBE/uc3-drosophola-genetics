import sys
from collections import defaultdict as d

from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--meta", dest="MA", help="Input file")

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

def countlist(lst):
    consec = [1]
    for x, y in zip(lst, lst[1:]):
        if x == y - 1:
            consec[-1] += 1
        else:
            consec.append(1)
    return consec

Header=["GapLength"]
for l in load_data(options.MA):
    if l.startswith("sample"):
        continue
    a=l.rstrip().split(",")
    Header.append(a[0])

print("\t".join(Header))

C=0
Cont=d(lambda: d(list))
for l in load_data(options.IN):
    a=l.rstrip().split()
    pops=a[3:]
    for i in range(len(pops)):
        if pops[i] == ".,.":
            Cont[a[0]][i].append(C)
    C+=1

Counts=d(lambda:d(int))
for Chrom,v in Cont.items():
    for pop,List in v.items():
        for i in countlist(List):
            Counts[i][pop]+=1
    
for Length,v in sorted(Counts.items()):
    PL=[]
    for i in range(len(pops)):
        PL.append(v[i]/C)
    print(str(Length)+"\t"+"\t".join([str(x) for x in PL]))
