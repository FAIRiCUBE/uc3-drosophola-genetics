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

geno_file = []
keep=[]
kick=[]
allfromvcf=[]
# columns=[*range(9,(n+9),1)]


for line in vcf:
    if line.startswith("#CHROM"):
        popcol = line.split()

with open(options.SAMPLES, 'r') as f:
    matching_pops = [line.strip() for line in f]
    columns=[i for i, x in enumerate(popcol) if x in matching_pops ]

meta = open(options.META, 'r').readlines()
popsize=[]
for line in meta:
    for pop in matching_pops:
        if line.startswith(pop):
            popsize.append(line.split(",")[4])

with open("results/BAYPASS/size.poolsize",'w') as file:
    file.write(' '.join(map(str, popsize)))


#columns = [9,11,13,15,16,17,18,21,23,25,29,31,32,33,36,37,41,42,44,47,49,50,56] # Selction only of populations we are going to use for spring
#columns = [10,12,14,19,20,22,28,34,35,38,39,40,46,48,51,52,53,54,55] # Selction only of populations we are going to use for autum
#columns = [*range(9,246+9,1)] # Selection of 246 populations

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
    for x in columns:
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
    sumref=0
    sumalt=0
    
    for j in range(0,len(geno_file_line),2):
        sumref=sumref+int(geno_file_line[j])
        #print(sumref)
    for k in range(1,len(geno_file_line),2):
        sumalt=sumalt+int(geno_file_line[k])
        #print(sumalt)
    if sumref!=0 and sumalt!=0:
        geno_file.append(geno_file_line)
        keep.append((chromosome+"\t"+position))
    else:
        kick.append((chromosome+"\t"+position))
    allfromvcf.append((chromosome+"\t"+position))



   



print("samples passed the monomorphic filter: " + str(len(geno_file)))
geno_output = open(options.OUT, 'w')

#for geno in geno_file:
#    print >> geno_output, "\t".join(geno)

for geno in geno_file:
    geno_output.write("\t".join(geno))
    geno_output.write("\n")


##write the poolsize file as well
#meta = open(options.META, "rt").readlines()

#write keep and kick

MarkersKept= open("results/BAYPASS/MarkersKept.csv", 'w')
for marker in keep:
    MarkersKept.write(marker)
    MarkersKept.write("\n")

MarkersKicked= open("results/BAYPASS/MarkersKicked.csv", 'w')
for marker in kick:
    MarkersKicked.write(marker)
    MarkersKicked.write("\n")

All= open("results/BAYPASS/AllMarkers.csv",'w')
for marker in allfromvcf:
    All.write(marker)
    All.write("\n")