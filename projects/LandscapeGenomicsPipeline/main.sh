## Author: Sonja Steindl
## Date: 23.02.2022
## Status: in progress 

## Load required Packages and Programs
## This Pipeline is set up to work with "modules" 
module load Tools/vcftools_0.1.13
module load Tools/BayPass_2.3

# Declare your work environment and provide a file with the names of the populations you want to analyze
wd=$(pwd)
cd $wd

#please provide a file with the names of the samples you want to include into the analysis and the path to it

#recoded samplenames due to new file structure 
#awk 'FNR==NR {file2[$1]; next } $2 in file2' data/samplenames.csv data/dest_v2.samps_25Feb2023.csv > data/samplenamesnew.csv
samples="data/samplenamesnew.csv"
#please set as input which chromosmal region you want to analyze


arm="3R"
mkdir results
mkdir data

# Workflow
# Starting Point Of This Pipeline: VCF-FORMAT 
# 1) Get (sample) dataset as VCF by:
cd data

#get the VCF
wget "http://berglandlab.uvadcos.io/vcf/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
mv dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz data/dest.PoolSeq.2023.ann.vcf.gz

#get wordclim data
wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/worldClim/dest.worldclim.csv
### get metadatafile
wget https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.csv

cd ..

# Declare the VCF File as Input File for the following steps
input="data/dest.PoolSeq.2023.ann.vcf.gz"
output="Subsampled_"${arm}".recode.vcf.gz" #name must match with the awk of the chromosomes 
outaf="Subsampled_"${arm}".recode.af"
scripts=${wd}"/scripts"


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
python3 ${scripts}/SubsampleVCF.py \
  --input results/${output} \
  --snps 10000 |
  gzip >results/k10.${output}

### convert to AFs
python3 ${scripts}/VCF2AF.py \
  --input results/k10.${output} \
  >results/k10.${outaf}

### restrict to samples and two biovariables7
## ATTENTION, column 24 added for nFLIES for baypass 
#this restriction only for glm?? later on py script 
#maybe make this "user friendlyier" and let user choose which variable to analyze

{
  head -1 data/dest_v2.samps_25Feb2023.csv
  grep -f ${samples} data/dest_v2.samps_25Feb2023.csv
} |
  cut -d "," -f1,4,5,6,7,24,30 \
    >results/metadata.csv

metadata="results/metadata.csv"

mkdir GM 
# Perform linear regression on 3R
#INCLUDE a Rmd File that captures statistics on the calculations and analysis 
Rscript ${scripts}/Plot_pvalues.R $wd results/k10.${outaf} $metadata 


# Latent Factor Mixed Model
LeaOut="results/LEA"
mkdir $LeaOut

# choose number of latent factors (nK) and number of i
# use subset of neutral snps to estimate nK
# ... maybe list of intronic snps 
nK=3
nR=3 # number of calculation repetitions for each factor
var="long"

for rep in $(seq 1 $nR)
do
Rscript ${scripts}/LEA_RunLFMM.R $LeaOut/${var}_run$rep ${wd}/results/k10.${outaf} ${wd}/${metadata} $var $nK 1 &
done

wait

Rscript ${scripts}/LEA_ZPcalc.R $LeaOut $nK $nR ${wd}/results/k10.${outaf} $var

# BAYPASS analyses
#parse inputcsv=${metadata}
#sh /media/inter/ssteindl/DEST/LanGen/BAYPASS/baypass_main.sh 
#Script "main" (Including geno_creation.py, some shell commands to create necessary files, run Baypass)
bayin="results/k10."${output}
baydir="results/BAYPASS"
mkdir $baydir
bayout=${baydir}"/baypass.geno"
baycov="results/BAYPASS/covariates.csv"

##note that this starts from VCF and not AF, needs to be changed in order to give comparable results 
###change in line 103 necessar<: output to out.af!!!
#python3 /media/inter/ssteindl/DEST/LanGen/BAYPASS/geno_creation_ext.py --input $input --output ${baydir}${bayout}
python3 ${scripts}/geno_creation_morph.py \
    --input $bayin \
    --output $bayout \
    --samples $samples \
    --metadata $metadata

## This script scans which columns (covariates) are numeric and only uses these for generating the BAYPASS covariate file and also creates and additional "covariate.info" file
## for indication which variables went into the analyses 
## IDEA: maybe recode variables that are string into numeric, also NAs not accepted yet ?
## change data/dest_v2... to the FINAL metadata/environemtnal data table  
python3 ${scripts}/create_cov.py --input data/dest_v2.samps_25Feb2023.csv --output $baycov --samples $samples

##Result
#create on results folder with subdirectories (results for each method) and one general "comparison result"??
g_baypass -npop $(wc -l $samples) \
-gfile $bayout \
-efile $baycov \
-outprefix ${baydir}/BayPass \
-poolsizefile $baydir/size.poolsize


echo """
### XtX statistics over the 10k SNPs
XtX=read.table("/media/inter/ssteindl/FC/test/BayPass_summary_pi_xtx.out",h=T)$M_XtX
pod.thresh=quantile(XtX,probs=0.95)
plotdf <- data.frame(x=c(1:10000), y=XtX)
PLOT <- ggplot(plotdf, aes(x=x, y=y)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=pod.thresh,col="blue",lty=2)+
  ylab("XtX") + theme_bw() + xlab("Genomic Position")

ggsave("XtXstats.png",
     PLOT,
     width=15,
     height=5)
""" > whatever.R

Rscript whatever.R

echo """
file <- read.table("/media/inter/ssteindl/FC/test/BayPass_summary_betai_reg.out",h=T)
file[file$COVARIABLE==1,]$eBPis -> Cov1P

PLOT <- ggplot(plotdf, aes(c(1:10000), y=Cov1P)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=pod.thresh,col="blue",lty=2)+
  ylab("eBPis (-log10P)") + theme_bw() + xlab("Genomic Position")

ggsave("eBPis.png",
       PLOT,
       width=15,
       height=5)
""" > cov.R

## idea:
# drawback to snps and overlap in snps via VENN diagram 
# test non-overlaps




