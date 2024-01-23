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
#from collections import defaultdict as d
#from optparse import OptionParser, OptionGroup

print(sys.argv[1])
vcf = gzip.open("dest.PoolSeq.2020.ann.vcf.gz", "rt", encoding="utf-8").readlines()[1:]
y="outpath"
geno_file = []
# columns=[*range(9,(n+9),1)]

#columns = [9,11,13,15,16,17,18,21,23,25,29,31,32,33,36,37,41,42,44,47,49,50,56] # Selction only of populations we are going to use for spring
#columns = [10,12,14,19,20,22,28,34,35,38,39,40,46,48,51,52,53,54,55] # Selction only of populations we are going to use for autum
columns = [*range(9,246+9,1)] # Selction of 246 populations

flag = 0
geno_output = open(y, 'w')

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
    geno_file.append(geno_file_line)
    flag = flag + 1
    print(flag)

print("longitude: " + str(len(geno_file)))
#geno_output = open("/media/inter/ssteindl/DEST/LanGen/BAYPASS/dummy/Biallelic.geno", 'w')

#for geno in geno_file:
#    print >> geno_output, "\t".join(geno)

for geno in geno_file:
    geno_output.write("\t".join(geno))
    geno_output.write("\n")
