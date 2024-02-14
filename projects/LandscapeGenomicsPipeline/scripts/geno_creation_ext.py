#!/usr/bin/python

# Script oroginated from: https://github.com/GonzalezLab
# Custom Changes by Sonja Steindl 
# 22.02.2023
# Creates geno file needed as input for running Baypass from VCF file


# Info: Multiallelic Sites will be reduced to biallelic, using the more frequent alternative allele
# Missing Data in the VCF (./.:.:.:.:.) is converted to "0" for the Baypass Program.
# Using 0 for both the REF and ALT allele therfore resembles missing data.
# This script is suitable to read .gz vcf files
# Number of populations in the vcf file need to be specified via "columns=[]"

import os
import sys
import csv
import subprocess
import gzip
import re
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="outfile")
parser.add_option("--samples", dest="SAMPLES", help="samplenames file")
parser.add_option("--metadata", dest="META", help="metadata file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

vcf = gzip.open(options.IN, "rt", encoding="utf-8").readlines()[1:]
#vcf = gzip.open("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/Subsample3.vcf.gz", "rt").readlines()[1:]
#vcf=gzip.open("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_4/results/2L/Subsampled_2L.recode.vcf.h5.gz", "rt", encoding="utf-8").readlines()[1:]
geno_file = []
#columns=[*range(9,(n+8),1)]

for line in vcf:
    if line.startswith("#CHROM"):
        popcol = line.split()

with open(options.SAMPLES, 'r') as f:
    matching_pops = [line.strip() for line in f]
    columns=[i for i, x in enumerate(popcol) if x in matching_pops ]
    print(matching_pops)
    #print(columns)

#with open("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/european.csv", 'r') as f:
#    matching_pops = [line.strip() for line in f]
#    columns=[i for i, x in enumerate(popcol) if x in matching_pops ]
#    print(matching_pops)
#    #print(columns)

meta = open(options.META, 'r').readlines()
header = meta[0]
meta_lines = meta[1:]
#meta= open("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/northamerica_meta.csv", 'r').readlines()
nFliesIndex=header.split(",").index("nFlies")
popsize=[]
for line in meta_lines:
    for pop in matching_pops:
        if line.startswith(pop):
            print(line)
    popsize.append(line.split(",")[nFliesIndex])

output_file_path = os.path.join(os.path.dirname(options.OUT), "size.poolsize")

# Write the data to the output file with error handling
try:
    with open(output_file_path, 'w') as file:
        file.write(' '.join(map(str, popsize)))
    print(f"File '{output_file_path}' created successfully.")
except Exception as e:
    print(f"Error writing to '{output_file_path}': {e}")


poolsize_double = os.path.join(os.path.dirname(options.OUT), "size_double.poolsize")

for i in range(len(popsize)):
    popsize[i] = int(popsize[i])*2

try:
    with open(poolsize_double, 'w') as file:
        file.write(' '.join(map(str, popsize)))
    print(f"File '{poolsize_double}' created successfully.")
except Exception as e:
    print(f"Error writing to '{poolsize_double}': {e}")


#with open("results/BAYPASS/size.poolsize",'w') as file:
#    file.write(' '.join(map(str, popsize)))


#columns = [9,11,13,15,16,17,18,21,23,25,29,31,32,33,36,37,41,42,44,47,49,50,56] # Selction only of populations we are going to use for spring
#columns = [10,12,14,19,20,22,28,34,35,38,39,40,46,48,51,52,53,54,55] # Selction only of populations we are going to use for autum
#columns = [*range(9,246+9,1)] # Selection of 246 populations
#columns=[0]
flag = 0
geno_output = open(options.OUT, 'w')

for line in vcf:
    if line.startswith("#"):
        continue
    geno_file_line = []
    chromosome = line.split("\t")[0]
    position = line.split("\t")[1]
    #geno_file_line.append(chromosome)
    #geno_file_line.append(position)
    for x in list(range(9,len(line.split("\t")))):
        population = line.split("\t")[x]
        alternative = population.split(":")[2]
        original = population.split(":")[1]
        #original = int(total) - int(alternative)
        if original != str("."):
            geno_file_line.append(str(original))
        else:
            geno_file_line.append(str("0"))
        if alternative != str("."):
            if bool(re.search(".*,.*", alternative)) == True: 
                a = alternative.split(",")[1]
                b = alternative.split(",")[0]
                geno_file_line.append(str(max(a,b)))
            else:
                geno_file_line.append(str(alternative))
        else: 
            geno_file_line.append(str("0"))
    geno_file.append(geno_file_line)
    flag = flag + 1
    #print(flag)

#print("SNP: " + str(len(geno_file)))
geno_output = open(options.OUT, 'w')

#for geno in geno_file:
#    print >> geno_output, "\t".join(geno)

for geno in geno_file:
    geno_output.write("\t".join(geno))
    geno_output.write("\n")

##write the poolsize file as well
#meta = open(options.META, "rt").readlines()


