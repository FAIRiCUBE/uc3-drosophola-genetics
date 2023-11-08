## Author: Sonja Steindl
## Status: in progress 02.10.2023


####### REQUIRED PROGRAMS #######
# VCFTOOLS   please install     #
# BAYPASS    please install     #
#################################

# Declare your work environment and provide a file with the names of the populations you want to analyze
wd=$(pwd)
cd $wd
mkdir data
mkdir results


module load Tools/vcftools_0.1.13

# Workflow
# Starting Point Of This Pipeline: VCF-FORMAT 
# Get (sample) dataset as VCF by:

# Please note, this is not the latest DEST data set
cd data
wget --tries=inf "http://berglandlab.uvadcos.io/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gdsdest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
wget "https://github.com/DEST-bio/DESTv2/blob/main/populationInfo/dest_v2.samps_25Feb2023.csv"
##take all samples
awk '{FS=","}{if (NR!=1) {print $1}}' samps_25Feb2023.csv > samplenames.csv
mv data/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz data/PoolSeq2023.vcf.gz


samples="data/samplenames.csv"

input="data/PoolSeq2023.vcf.gz"
samplefile="data/samps_25Feb2023.csv"
### restrict analysis to samples and two biovariables7
    ## ATTENTION, colum 12 added for nFLIES for baypass
{
  head -1 $samplefile
  grep -f $samples $samplefile
} |
  cut -d "," -f1,4,5,12,6,7,24,30,41 \
    >results/metadata.csv

metadata="results/metadata.csv"

#genomic_regions=($(gunzip -c data/dest.PoolSeq.2023.norep.vcf.gz | awk '!/^#/ {print $1}' | uniq))
genomic_regions=($(gunzip -c data/PoolSeq2023.vcf.gz | awk '/^[^#]/ {print $1}' | sort -u))


for value in "${genomic_regions[@]}"; do
    arm=$value
    mkdir results/${arm}
    output="results/${arm}/Subsampled_"${arm}".recode.vcf.gz" #name must match with the awk of the chromosomes 
    outaf="results/${arm}/Subsampled_"${arm}".recode.af"

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
      gzip > ${output}

    ### randomly pick 10k lines
    python scripts/SubsampleVCF.py \
      --input ${output} \
      --snps 10000 |
      gzip >${output}.k10.gz

    ### convert to AFs
    python scripts/VCF2AF.py \
      --input ${output}.k10.gz \
      >${outaf}.k10

    # Perform linear regression on 3R
    Rscript scripts/Plot_pvalues.R $wd ${outaf}.k10 $metadata $arm
    
    
    # BAYPASS analyses
    #parse inputcsv=${metadata}
    #sh /media/inter/ssteindl/DEST/LanGen/BAYPASS/baypass_main.sh 
    #Script "main" (Including geno_creation.py, some shell commands to create necessary files, run Baypass)
    bayin="${output}.k10.gz"
    baydir="results/${arm}/BAYPASS/"
    mkdir $baydir
    bayout=${baydir}"baypass.geno"
    baycov=${baydir}"covariates.csv"

    ##note that this starts from VCF and not AF, needs to be changed in order to give comparable results 
    ###change in line 103 necessary: output to outaf!!!
    #python3 /media/inter/ssteindl/DEST/LanGen/BAYPASS/geno_creation_ext.py --input $input --output ${baydir}${bayout}
    python3 scripts/geno_creation_ext.py \
        --input $bayin\
        --output $bayout \
        --samples $samples \
        --metadata $metadata

        ##inlcude: if sum of alt/ref alleles = 0; continue
        ##and write file with pos of (non-continue blabla)

    ## IMPORTANT, MALE is a covariate which cannot be interpreted by BAYPASS (String?) and therefore recode?
    ## 
    python3 scripts/create_cov.py --input $samplefile --output $baycov --samples $samples
    genofile=${bayout}
    variables=$(<"results/2L/BAYPASS/covariates.covariate.info.csv")

  # Use a for loop to iterate over words (assuming space-separated words)
    #genofile="/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/POLYmorphic_baypass.geno"
    ##Result
    #create on results folder with subdirectories (results for each method) and one general "comparison result"??
    /media/inter/ssteindl/DEST/LanGen/baypass_2.3/sources/g_baypass -npop $(wc -l $samples) -gfile $genofile -efile $baycov -outprefix $baydir/BayPass -poolsizefile ${baydir}/size.poolsize

    Rscript scripts/baypass_plotXtX.R ${baydir}/BayPass_summary_pi_xtx.out ${arm}

    Rscript scripts/baypass_plot_eBPis.R ${baydir}/BayPass_summary_betai_reg.out ${arm} results/2L/BAYPASS/covariates.covariate.info.csv
  

    # Latent Factor Mixed Model
    LeaOut="results/${arm}/LEA"
    mkdir $LeaOut

    # Choose number of estimated latent factors (nK) and number of i
    # Number of calculation repetitions for each factor.
    nR=3
    nK=7

    for cov in $variables
    do
    ## IMPORTANT: make the variable iterable 
      var=$cov
      for rep in $(seq 1 $nR)
      do
      echo $wd
      echo $metadata
      echo $rep
      Rscript scripts/LEA_RunLFMM.R $LeaOut/${var}_run$rep ${wd}/${outaf}.k10 ${wd}/${metadata} $var 1 
      done
    Rscript scripts/LEA_ZPcalc.R $LeaOut $nK $nR ${wd}/${outaf}.k10 $var
  done
done





