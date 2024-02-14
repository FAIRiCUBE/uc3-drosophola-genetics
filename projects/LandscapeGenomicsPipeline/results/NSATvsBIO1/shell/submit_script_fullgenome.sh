#!/bin/sh
## name of Job
#PBS -N LandscapeGenomicsfullgenome
## Redirect output stream to this file.
#PBS -o /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/logs
## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe
## Select a maximum of 20 cores and 200gb of RAM
#PBS -l select=1:ncpus=20:mem=200gb
## load all necessary software into environment

cd /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9

module load Tools/vcftools_0.1.13
     
mkdir results/fullgenome
output="results/fullgenome/Subsampled_fullgenome.recode.vcf" #name must match with the awk of the chromosomes 
outaf="results/fullgenome/Subsampled_fullgenome.af"

### This step removes polyploidies, focus on region (Chromosome) and subsample 9 population samples and exlcude all sites with missing data
zcat /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/PoolSeq2023.vcf.gz |
  awk '$0~/^\#/ || length($5)==1' |
  vcftools \
    --vcf - \
    --keep /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/data/samplelist.csv \
    --stdout \
    --recode |
  grep -v "\./\." |
  gzip > ${output}.gz

### randomly pick 10k lines
python /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/SubsampleVCF.py \
  --input ${output}.gz \
  --snps 10000 |
  gzip >${output}.h5.gz

### convert to AFs
python /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/VCF2AF.py \
  --input ${output}.h5.gz \
  >${outaf}.h5

# Perform linear regression on 3R
Rscript ../scripts/Plot_pvalues.R /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9 ${outaf}.h5 data/metadata.csv fullgenome

# BAYPASS analyses

#Script "main" (Including geno_creation.py, some shell commands to create necessary files, run Baypass)
bayin="${output}.h5.gz"
baydir="results/fullgenome/BAYPASS"
mkdir $baydir
bayout=${baydir}"/baypass.geno"
baycov=${baydir}"/covariates.csv"

python3 /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/geno_creation_ext.py     --input $bayin \
    --output $bayout \
    --samples /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/data/samplelist.csv \
    --metadata data/metadata.csv
    ##inlcude: if sum of alt/ref alleles = 0; continue
    ##and write file with pos of (non-continue blabla)

## IMPORTANT, MALE is a covariate which cannot be interpreted by BAYPASS (String?) and therefore recode?

python3 /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/create_cov.py --input data/metadata.csv --output ${baycov} --samples /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/data/samplelist.csv
genofile=${bayout}
variables=$(<"results/fullgenome/BAYPASS/covariates.covariate.info.csv")
  #Use a for loop to iterate over words (assuming space-separated words)

/media/inter/ssteindl/DEST/LanGen/baypass_2.3/sources/g_baypass -npop 3 -gfile $genofile -efile $baycov -outprefix $baydir/BayPass -poolsizefile ${baydir}/size.poolsize
Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/baypass_plotXtX.R ${baydir}/BayPass_summary_pi_xtx.out fullgenome
Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/baypass_plot_eBPis.R ${baydir}/BayPass_summary_betai_reg.out fullgenome results/fullgenome/BAYPASS/covariates.covariate.info.csv

# Latent Factor Mixed Model
LeaOut="results/fullgenome/LEA"
mkdir $LeaOut

# Choose number of estimated latent factors (nK) and number of i
# Number of calculation repetitions for each factor.
nR=3
nK=7

#outaf="results/fullgenome/Subsampled_fullgenome.af"

IFS=' ' read -ra variables < "results/fullgenome/BAYPASS/covariates.covariate.info.csv"

#for variable in "${variables[@]}"
#do
#var=$variable
#  for rep in $(seq 1 $nR)
#  do
#  echo /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9
#  echo data/metadata.csv
#  echo $rep
#  echo $var
#  Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/LEA_RunLFMM.R ${wd}/${LeaOut}/${var}_run$rep /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/$outaf.h5 /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/data/metadata.csv $var 1 
#  done
#  Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/LEA_ZPcalc.R $LeaOut $nK $nR /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/$outaf.h5 $var
#done

for variable in "${variables[@]}"
do
var=$variable
  for rep in $(seq 1 $nR)
  do
  echo /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9
  echo data/metadata.csv
  echo $rep
  echo $var
  Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/LEA_RunLFMM.R $wd/$LeaOut/${var}_run$rep /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/$outaf.h5 /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/data/metadata.csv $var 1 
  done
  Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/LEA_ZPcalc.R $wd/$LeaOut $nK $nR /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_9/$outaf.h5 $var
done


##from all adjusted p values filter and intersect


