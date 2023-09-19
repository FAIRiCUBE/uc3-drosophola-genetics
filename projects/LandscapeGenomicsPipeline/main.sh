## Author: Sonja Steindl
## Status: in progress 


####### REQUIRED PROGRAMS #######
# VCFTOOLS   please install     #
# BAYPASS    please install     #
#################################

# Declare your work environment and provide a file with the names of the populations you want to analyze
wd=$(pwd)
cd $wd
#mkdir data
#please provide a file with the names of the samples you want to include into the analysis and the path to it
samples="samplenames.csv"
#samples="/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/data/samps_10Nov2020.csv"
#please set as input which chromosmal region you want to analyze
arm="3R"
mkdir results
cd results 
mkdir $arm
cd $wd

## Load required Packages and Programs
## This Pipeline is set up to work with "modules" 
module load Tools/vcftools_0.1.13

# Workflow
# Starting Point Of This Pipeline: VCF-FORMAT 
# Get (sample) dataset as VCF by:
# Please note, this is not the latest DEST data set 
wget "http://berglandlab.uvadcos.io/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
### get metadatafile
wget https://raw.githubusercontent.com/DEST-bio/DEST_freeze1/main/populationInfo/samps_10Nov2020.csv
mv dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz data/dest.PoolSeq.2020.ann.vcf.gz


# Declare the VCF File as Input File for the following steps
input="data/dest.PoolSeq.2020.ann.vcf.gz"
output="results/${arm}/Subsampled_"${arm}".recode.vcf.gz" #name must match with the awk of the chromosomes 
outaf="results/${arm}/Subsampled_"${arm}".recode.af"
#scripts="/media/inter/mkapun/projects/FAIRiCUBE/uc3-drosophola-genetics/DrosoEnvironmentTest/scripts"


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



### restrict to samples and two biovariables7
## ATTENTION, colum 12 added for nFLIES for baypass
{
  head -1 samps_10Nov2020.csv
  grep -f ${samples} samps_10Nov2020.csv
} |
  cut -d "," -f1,5,12,6,7,30,41 \
    >results/metadata.csv

metadata="results/metadata.csv"

# Perform linear regression on 3R
Rscript /media/inter/ssteindl/DEST/LanGen/Plot_pvalues.R $wd results/${outaf} $metadata
# Latent Factor Mixed Model
LeaOut="results/${arm}/LEA"
mkdir $LeaOut
# choose number of estimated latent factors (nK) and number of i
nK=3
#number of calculation repetitions for each factor
nR=3 

## IMPORTANT: make the variable iterable 
var="lat"
for rep in $(seq 1 $nR)
do
Rscript /media/inter/ssteindl/DEST/LanGen/LEA_RunLFMM.R $LeaOut/${var}_run$rep ${wd}/results/${outaf} ${wd}/${metadata} $var $nK 1 &
done
wait

Rscript /media/inter/ssteindl/DEST/LanGen/LEA_ZPcalc.R $LeaOut $nK $nR ${wd}/results/${outaf} $var

# BAYPASS analyses
#parse inputcsv=${metadata}
#sh /media/inter/ssteindl/DEST/LanGen/BAYPASS/baypass_main.sh 
#Script "main" (Including geno_creation.py, some shell commands to create necessary files, run Baypass)
bayin="results/k10."${output}
baydir="results/BAYPASS/"
mkdir $baydir
bayout=${baydir}"baypass.geno"
baycov=${baydir}"covariates.csv"

##note that this starts from VCF and not AF, needs to be changed in order to give comparable results 
###change in line 103 necessary: output to outaf!!!
#python3 /media/inter/ssteindl/DEST/LanGen/BAYPASS/geno_creation_ext.py --input $input --output ${baydir}${bayout}
python3 /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/geno_creation_polymorphic.py \
    --input /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/k10.Subsampled_3R.recode.vcf.gz\
    --output $bayout \
    --samples $samples \
    --metadata $metadata

    ##inlcude: if sum of alt/ref alleles = 0; continue
    ##and write file with pos of (non-continue blabla)

## IMPORTANT, MALE is a covariate which cannot be interpreted by BAYPASS (String?) and therefore recode?
## 
python3 /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/create_cov.py --input samps_10Nov2020.csv --output $baycov #--samples $samples

genofile=${bayout}

genofile="/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/POLYmorphic_baypass.geno"
##Result
#create on results folder with subdirectories (results for each method) and one general "comparison result"??
/media/inter/ssteindl/DEST/LanGen/baypass_2.3/sources/g_baypass -npop $(wc -l $samples) -gfile $genofile -efile $baycov -outprefix results/BAYPASS/Polymorph/BayPass -poolsizefile results/BAYPASS/size.poolsize

echo """
### XtX statistics over the 10k SNPs
library(ggplot2)
XtX=read.table(\"/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/BayPass_summary_pi_xtx.out\",h=T)$M_XtX
pod.thresh=quantile(XtX$XtXst,probs=0.95)
plotdf <- data.frame(x=c(1:10000), y=XtX)
PLOT <- ggplot(plotdf, aes(x=x, y=y.XtXst)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=pod.thresh,col=\"blue\",lty=2)+
  ylab(XtX) + theme_bw() + xlab(\"Genomic Position\")
ggsave(\"XtXstats.png\",
     PLOT,
     width=15,
     height=5)
""" > /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/intermediate/Stats.R

Rscript /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/intermediate/Stats.R

echo """
file <- read.table(\"/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/BayPass_summary_betai_reg.out\",h=T)
file[file$COVARIABLE==1,]$eBPis -> Cov1P

PLOT <- ggplot(plotdf, aes(c(1:10000), y=Cov1P)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=pod.thresh,col=\"blue\",lty=2)+
  ylab(\"eBPis (-log10P)\") + theme_bw() + xlab(\"Genomic Position\")

ggsave(\"eBPis.png\",
       PLOT,
       width=15,
       height=5)
""" > /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/intermediate/cov.R

