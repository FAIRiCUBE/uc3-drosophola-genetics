#/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/PoolSeq2023.vcf.gz

## Author: Sonja Steindl
## Status: in progress 02.10.2023


####### REQUIRED PROGRAMS #######
# VCFTOOLS   please install     #
# BAYPASS    please install     #
#################################



DIRECTORY=$(pwd)

while getopts d:o:s:i:m:w:v:b:? opt 2>/dev/null
do
  case $opt in
    d) DIRECTORY=$OPTARG;;
    o) OUTPUTDIR=$OPTARG;;
    s) SAMPLES=$OPTARG;;
    i) VCFFILE=$OPTARG;;
    m) METADATA=$OPTARG;;
    w) WORLDCLIM=$OPTARG;;
    v) VARIABLE=$OPTARG;;
    b) BIOVAR=$OPTARG;;
    ?) echo "Valid parameters are: [-l, -s]" ;;
  esac
done

if [ -z "$DIRECTORY" ] || [ -z "$SAMPLES" ] || [ -z "$VCFFILE" ]; then
  echo "Error: Missing required parameters. Please provide values for -d, -s, and -i."
  exit 1
fi

if [ -n "$DIRECTORY" ]; then
  echo "DIRECTORY is set to: $DIRECTORY"
fi

if [ -n "$SAMPLES" ]; then
  echo "SAMPLES is set to: $SAMPLES"
fi

if [ -n "$VCFFILE" ]; then
  echo "VCFFILE is set to: $VCFFILE"
fi

if [ -n "$METADATA" ]; then
  echo "METADATA is set to: $METADATA"
fi

if [ -n "$WORLDCLIM" ]; then
  if [ -n "$BIOVAR" ]; then
    echo "METADATA is set to: $METADATA. Including Parameter BIOVAR: $BIOVAR"
  else
    echo "Error: BIOVAR is missing. Please set the BIOVAR parameter."
    exit 1
  fi
fi

if [ -n "$VARIABLE" ]; then
  echo "VARIABLE is set to: $VARIABLE"
fi

scriptdir=$(pwd)/scripts
variable=$VARIABLE
metadata_dest=$METADATA
samples=$SAMPLES
input=$VCFFILE
#wd=$(pwd)
wd=$DIRECTORY/$OUTPUTDIR
mkdir $wd
cd $wd

mkdir data
mkdir results

### restrict analysis to samples and two biovariables7
    ## ATTENTION, colum 12 added for nFLIES for baypass

#npop=$(wc -l < "$samples")
echo "PREPARING METADATA FILE"

python3 ${scriptdir}/MergeData.py \
  --metadata ${metadata_dest} \
  --variable ${variable} \
  --samplenames ${samples} \
  --worldclim $WORLDCLIM \
  --biovariable $BIOVAR \
  --output data/metadata.csv

echo "MERGING DONE"

metadata=data/metadata.csv

npop=$(tail -n +2 ${metadata} | wc -l)

###overwrite samples to what is available in data
#samples=$(awk -F ',' 'NR>1 {print $1}' data/metadata.csv)
awk -F ',' 'NR>1 {print $1}' $metadata > data/samplelist.csv

samplelist="${wd}/data/samplelist.csv"


echo "READING GENOMIC REGIONS DONE"

shellfolder="${wd}/shell"
logfolder="${wd}/logs"
mkdir $logfolder
mkdir $shellfolder


arm="fullgenome"
echo $arm
script_file="submit_script_${arm}.sh"
cat > "${shellfolder}/$script_file" <<EOL
#!/bin/sh
## name of Job
#PBS -N LandscapeGenomics${arm}
## Redirect output stream to this file.
#PBS -o ${logfolder}
## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe
## Select a maximum of 20 cores and 200gb of RAM
#PBS -l select=1:ncpus=20:mem=200gb
## load all necessary software into environment

cd ${wd}

module load Tools/vcftools_0.1.13
    
mkdir results/${arm}
output="results/${arm}/Subsampled_${arm}.recode.vcf" #name must match with the awk of the chromosomes 
outaf="results/${arm}/Subsampled_${arm}.af"

### This step removes polyploidies, focus on region (Chromosome) and subsample 9 population samples and exlcude all sites with missing data
zcat ${input} |
  awk '\$0~/^\#/ || length(\$5)==1' |
  vcftools \\
    --vcf - \\
    --keep ${samplelist} \\
    --stdout \\
    --recode |
  grep -v "\./\." |
  gzip > \${output}.gz

### randomly pick 10k lines
python ${scriptdir}/SubsampleVCF.py \\
  --input \${output}.gz \\
  --snps 10000 |
  gzip >\${output}.h5.gz

### convert to AFs
python ${scriptdir}/VCF2AF.py \\
  --input \${output}.h5.gz \\
  >\${outaf}.h5

# Perform linear regression on 3R
Rscript ../scripts/Plot_pvalues.R $wd \${outaf}.h5 $metadata $arm

# BAYPASS analyses

#Script "main" (Including geno_creation.py, some shell commands to create necessary files, run Baypass)
bayin="\${output}.h5.gz"
baydir="results/${arm}/BAYPASS"
mkdir \$baydir
bayout=\${baydir}"/baypass.geno"
baycov=\${baydir}"/covariates.csv"

python3 ${scriptdir}/geno_creation_ext.py \
    --input \$bayin \\
    --output \$bayout \\
    --samples $samplelist \\
    --metadata $metadata
    ##inlcude: if sum of alt/ref alleles = 0; continue
    ##and write file with pos of (non-continue blabla)

## IMPORTANT, MALE is a covariate which cannot be interpreted by BAYPASS (String?) and therefore recode?

python3 /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/create_cov.py --input ${metadata} --output \${baycov} --samples $samplelist
genofile=\${bayout}
variables=\$(<"results/${arm}/BAYPASS/covariates.covariate.info.csv")
  #Use a for loop to iterate over words (assuming space-separated words)

/media/inter/ssteindl/DEST/LanGen/baypass_2.3/sources/g_baypass -npop $npop -gfile \$genofile -efile \$baycov -outprefix \$baydir/BayPass -poolsizefile \${baydir}/size.poolsize
Rscript ${scriptdir}/baypass_plotXtX.R \${baydir}/BayPass_summary_pi_xtx.out ${arm}
Rscript ${scriptdir}/baypass_plot_eBPis.R \${baydir}/BayPass_summary_betai_reg.out ${arm} results/${arm}/BAYPASS/covariates.covariate.info.csv

# Latent Factor Mixed Model
LeaOut="results/${arm}/LEA"
mkdir \$LeaOut

# Choose number of estimated latent factors (nK) and number of i
# Number of calculation repetitions for each factor.
nR=3
nK=7
for cov in \$variables
do
## IMPORTANT: make the variable iterable 
  var=\$cov
  for rep in \$(seq 1 \$nR)
  do
  echo $wd
  echo $metadata
  echo \$rep
  Rscript ${scriptdir}/LEA_RunLFMM.R \${wd}/\${LeaOut}/\${var}_run\$rep $wd/\$outaf.h5 ${wd}/${metadata} \$var 1 
  done
Rscript ${scriptdir}/LEA_ZPcalc.R \$LeaOut \$nK \$nR ${wd}/\$outaf.h5 \$var
  done
EOL
qsub "${shellfolder}/$script_file"
echo "DONE"





