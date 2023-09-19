# Landscape Genomics Pipeline 

This repository is supposed to hold documentation on the DEST_Pipeline used for the UC3 Project of FAIRiCUBE. 

## Objectives

This pipeline performs several landscape-genomics analysis on a SNP called data set. 

## Preliminary Results on Test Data


## Requirements
- VCFTools (vcftools.github.io/license.html)
- Baypass2.3 (http://www1.montpellier.inra.fr/CBGP/software/baypass/index.html)
- LEA R-Package (Citation)

## Workflow Contents 

0) [Set Up Environment](#set-up-environment)
1) [Download Data](#download-data)
2) [Subsample VCF for chromosomal arms](#subsample-vcf-for-chromosome)

2.1) [Remove Polyploides](#subsample-vcf-for-chromosome)

2.2) [Subsample 10k Variants](#subsample-vcf-for-chromosome)

2.3) [Convert to Allele Frequencies](#subsample-vcf-for-chromosome)

3) [Define Metadata / Covariates for analyses](#define-metadata) 

4)  [Linear Model](#linear-model)

5) [Latent Factor Mixed Model](#latent-factor-mixed-model)

6) [Baypass Analysis](#baypass-analysis)

**Creation of geno file**: Keep in mind, for this pipeline only a few samples are used therefore there are more likely "monomorphic" sites, if there are no monomorphic sites all will be kept and all positions are represented in the keep-file.



7) [Result Comparison](#Result-Comparison)

x) The basic workflow can be performed according to [main.sh](d/d/main.sh)

---

## Workflow Step by Step

### Set Up Environment

```bash
# Declare your work environment and provide a file with the names of the populations you want to analyze
cd $wd
mkdir data

#provide a file with the names of the samples you want to include into the analysis 
samples="samplenames.csv"

#set as input which chromosmal region you want to analyze
arm="3R"
mkdir results
mkdir results/${arm}
cd $wd

#Load required Packages and Programs
module load Tools/vcftools_0.1.13
```

### Download Data

```bash
#Download via wget
wget "http://berglandlab.uvadcos.io/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
wget https://raw.githubusercontent.com/DEST-bio/DEST_freeze1/main/populationInfo/samps_10Nov2020.csv
mv dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz data/dest.PoolSeq.2020.ann.vcf.gz
mv samps_10Nov2020.csv data/samps_10Nov2020.csv

# Declare the VCF File as Input File for the following steps
input="data/dest.PoolSeq.2020.ann.vcf.gz"
output="results/${arm}/Subsampled_"${arm}".recode.vcf.gz"
outaf="results/${arm}/Subsampled_"${arm}".recode.af"

```


### Subsample VCF for chromosome

```bash
### This step removes polyploidies, focus on 3R (Chromosome) and subsample 9 population samples and exlcude all sites with missing data
zcat ${input} |
  awk '$0~/^\#/ || length($5)==1' |
  awk '$0~/^\#/ || $1=="'$arm'" ' |
  vcftools \
    --vcf - \
    --keep ${samples} \
    --stdout \
    --recode |
  grep -v "\./\." |
  gzip >results/${output}

### randomly pick 10k lines
python ${scripts}/SubsampleVCF.py \
  --input results/${output} \
  --snps 10000 |
  gzip >results/k10.${output}

### convert to AFs
python ${scripts}/VCF2AF.py \
  --input results/k10.${output} \
  >results/k10.${outaf}

```

### Define Metadata

```bash
### restrict to samples and two biovariables7
## ATTENTION, colum 12 added for nFLIES for baypass
{
  head -1 samps_10Nov2020.csv
  grep -f ${samples} samps_10Nov2020.csv
} |
  cut -d "," -f1,5,12,6,7,30,41 \
    >results/metadata.csv

metadata="results/metadata.csv"
```