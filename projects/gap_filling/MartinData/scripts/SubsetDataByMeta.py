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

## Define criteria and arguments for filtering
Criteria=options.cri.split("=")[0]
Arguments = options.cri.split("=")[1].split(",")

C=0
Keep=[]
## open META file
O1=open(options.out+".meta","wt")
for l in load_data(options.meta):
    if l.startswith("sample"):
        ## get header
        header=l.rstrip().split(",")
        continue
    ## Make dictionary based on header for every entry
    DATA = dict(zip(header, l.rstrip().split(",")))

    ## only retain entries, i.e., Row IDs, if criteria == arguments
    if DATA[Criteria] in Arguments:
        Keep.append(C)
        O1.write(l)
    C+=1

O2 = gzip.open(options.out+".af.gz", "wt")
for l in load_data(options.IN):
    a=l.rstrip().split()

    ## keep Chrom, Pos and Alleles
    PL=a[:3]

    ## get columns with allele Counts
    pops=a[3:]

    ## only keep columns to keep based on criteria matching arguments in metadata
    KL=[pops[x] for x in Keep]

    ## sum up the counts for each allele across all pops- YEAH, list comprehensions
    TEST=[sum([int(z) for z in x if z!="."]) for x in zip(*[y.split(",") for y in KL])]

    ## test if SNP is polymorphic across all kept populations, i.e. all alleles with counts > 0, otherwise continue
    if 0 in TEST:
        continue

    ## write alleles counts to new file
    O2.write("\t".join(PL)+"\t"+"\t".join(KL)+"\n")
