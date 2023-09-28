# Landscape Genomics Pipeline 

This repository is supposed to hold documentation on the DEST_Pipeline (DROSEU) derivate used for the UC3 Project of FAIRiCUBE.
This is a standalone repository, as well as a directory within the FAIRiCUBE uc3-drosophila-genetics.
All the scripts used in the workflow are provided via the scripts directory.

## Objectives

This pipeline performs several landscape-genomics analysis on a SNP called data set intersected with environemntal data.
In this repository, the pipeline is documented for execution on 9 test populations, to easy computational workload in order to get familiar with the process.

![Flowchart](LandscapeGenomicsFairiCube.drawio.png)


## Requirements
- VCFTools
  - http://vcftools.github.io/license.html

- Baypass2.3 
  - http://www1.montpellier.inra.fr/CBGP/software/baypass/index.html

- LEA R-Package (Citation)
  - https://www.bioconductor.org/packages/release/bioc/html/LEA.html

    
- Working with Modules
  -    https://modules.readthedocs.io/en/stable/INSTALL.html


## Workflow Contents 

<!--ts-->
* [Set Up Environment](#set-up-environment)
* [Download Data](#download-genomic-data)
  * [Get Environmental Data](#download-environmental-data) 
* [Define Metadata / Covariates for analyses](#define-metadata-and-chromsomal-regions) 
* [Subsample VCF for chromosomal arms](#subsample-vcf-for-chromosome)
  * [Remove Polyploides](#subsample-vcf-for-chromosome)
  * [Subsample 10k Variants](#subsample-vcf-for-chromosome)
  * [Convert to Allele Frequencies](#subsample-vcf-for-chromosome)
* [Linear Model](#linear-model-with-r)
* [Latent Factor Mixed Model](#latent-factor-mixed-model)
* [Baypass Analysis](#baypass-analysis)
* [Result Interpretation and Comparison](#result-interpretation-and-comparison)
<!--te-->

x) The detailed workflow can be viewed and executed according to [main.sh](main.sh)


## Before Exploring The Pipeline 

Keep in mind, that the concept of the pipeline is to work with big data sets and therefore requires a lot of computational resources. We therfore provide a documentation and guideline to run the pipeline on a subset of population data of Drosophila melanogaster, provided by DEST (Drosophila Evolution over Space and Time) Data from DROSEU. Additionaly, structured guidance on how to run the pipeline with HPC scheduling systems will be included in this documentation.



## Workflow Step by Step

### Set Up Environment

This section does not contain specific guidance on installation of required programs. The purpose of this chapter is to define a working directory for the pipeline to run without any issues and to produce a structured output.
More information on how to install required software can be found in the [Requirements section](#requirements).

```bash
# Declare your work environment and provide a file with the names of the populations you want to analyze.
wd=$(pwd)
cd $wd
mkdir data
mkdir results
#Load required Packages and Programs.
module load Tools/vcftools_0.1.13
```

### Download Genomic Data

The data on geolocation of Drosohila melanogaster samples togehter with the according metadata, as well as the genetic data as VCF can be retrieved via DEST.bio. 

```bash
#Download DEST Drosophila data via wget as well as metadata created by DEST worling group.
cd data
wget --tries=inf "http://berglandlab.uvadcos.io/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
wget "https://raw.githubusercontent.com/DEST-bio/DEST_freeze1/main/populationInfo/samps_10Nov2020.csv"
mv dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz dest.PoolSeq.2020.ann.vcf.gz
samples="data/samps_10Nov2020.csv"
cd $wd
```

### Download Environmental Data

If you want to work with customly aquired environmental data, make sure you put it in the corresponding format, the number of variables (columns) is not limited. The file structure can be taken from [metadata.csv](results/metadata.csv) which will be created in the next step of the workflow. 

| Sample ID       | Variable 1      | Variable 2      |
| --------------- | --------------- | --------------- |
| AT_Kar_See_1_2014-08-17 | 44 | 255 |
| CH_Vau_Vul_1_2020-09-15| -128| 240 |
|CM_Nor_Oku_1_2004-04-15|  63|  54|
|  |  |  |




#### Define Metadata And Chromsomal Regions

```bash
#Create metadata.csv from the original DEST sample.csv
{
  grep -f ${samples} data/samps_10Nov2020.csv
} |
  cut -d "," -f1,5,12,6,7,30,41 \
    >results/metadata.csv

#Use predefined metadata or enter path to your custom metadata file
metadata="results/metadata_v2.csv"
```

Be aware, that ALL the following steps are all performed in a for-loop that is absed on chromosomal regions. The Loop continues after this section (see main.sh).

```bash
#Get all chromosomal regions in your data and set as input what you want to analyze. In the case of the provided VCF file by DEST.bio, regions are provided as chromosomal arms and analysis will be performed according to these arms.

genomic_regions=($(gunzip -c data/dest.PoolSeq.2023.norep.vcf.gz | awk '!/^#/ {print $1}' | uniq))
```

```bash
for value in "${genomic_regions[@]}"; do
    arm=$value
    mkdir results/${arm}
  
    # Declare the VCF File as Input File for the    following steps:
    input="data/dest.PoolSeq.2023.norep.vcf.gz"
    output="results/${arm}/Subsampled_"${arm}".recode.vcf.gz"
    outaf="results/${arm}/Subsampled_"${arm}".recode.af"

    #Provide a file with the names of the samples   you   want to include into the analysis, otherwise  use all  samples provided as default.
    samples="data/samps_10Nov2020.csv"
    #...

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

### Linear Model With R 

This is a simple regression statistics performed in R.

```bash
    Rscript /media/inter/ssteindl/DEST/LanGen/Plot_pvalues.R $wd results/${outaf} $metadata
```

### Latent Factor Mixed Model 

For this part of the analysis, two R scripts are performing the lfmm analysis and intermediate files (genotypes.lfmm & gradients.lfmm) are created.
Each iteration for analysing a variable creates a new folder.

```bash
    # Latent Factor Mixed Model

    LeaOut="results/${arm}/LEA"
    mkdir $LeaOut

    # Number of estimated latent factors (nK) and number of i is calculated within LEA_RunLFMM.R with the quick.elbow() function of bigpca-package.
    # To run quick.elbow(), manually install: https://cran.r-project.org/src/contrib/Archive/bigpca/bigpca_1.1.tar.gz

    # Number of calculation repetitions for each factor.
    nR=3 

    ## IMPORTANT: make the variable iterable 
    var="lat"
    for rep in $(seq 1 $nR)
    do
    Rscript scripts/LEA_RunLFMM.R $LeaOut/${var}_run1 ${wd}/results/${outaf} ${wd}/${metadata} "lat" 1 &
    done
    
    wait

    Rscript /media/inter/ssteindl/DEST/LanGen/LEA_ZPcalc.R $LeaOut $nK $nR ${wd}/results/${outaf} $var
```


### Baypass Analysis

**`Attention when creating geno file:`**  Keep in mind, for this pipeline only a few samples are used therefore there are more likely "monomorphic" sites, if there are no monomorphic sites all will be kept and all positions are represented in the keep-file.
In this analysis, a .geno and a .poolsize file are needed and created. 

```bash
    bayin="results/k10."${output}
    baydir="results/BAYPASS/"
    mkdir $baydir
    bayout=${baydir}"baypass.geno"
    baycov=${baydir}"covariates.csv"

    #python3 /media/inter/ssteindl/DEST/LanGen/BAYPASS/geno_creation_ext.py --input $input --output ${baydir}$    {bayout}
    python3 /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/geno_creation_polymorphic.py \
        --input /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/k10.Subsampled_3R.recode.vcf.gz\
        --output $bayout \
        --samples $samples \
        --metadata $metadata

    python3 /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/create_cov.py --input samps_10Nov2020.   csv --output $baycov #--samples $samples

    genofile=${bayout}

    genofile="/media/inter/ssteindl/FC/LandscapeGenomicsPipeline/results/BAYPASS/POLYmorphic_baypass.geno"
    ##Result
    #create on results folder with subdirectories (results for each method) and one general "comparison     result"??
    /media/inter/ssteindl/DEST/LanGen/baypass_2.3/sources/g_baypass -npop $(wc -l $samples) -gfile $genofile    -efile $baycov -outprefix results/BAYPASS/Polymorph/BayPass -poolsizefile results/BAYPASS/size.poolsize

    Rscript /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/baypass_plotXtX.R results/BAYPASS/   BayPass_summary_pi_xtx.out

    Rscript /media/inter/ssteindl/FC/LandscapeGenomicsPipeline/scripts/baypass_plot_eBPis.R results/BAYPASS/    BayPass_summary_betai_reg.out

done

```

### Result Interpretation And Comparison

As output of these analysis, we mainly observe -log10 p-values for each position in the genomic region. In a Manhattan Plot, genomic positions are plotted against -log10 p-values (mostly with ggplot). 


#### GM

In the worked example, linear regression for two biovariables was performed in different chromsomal arms. 
Below the results for loniguted and min_month performed on 2L are depicted.

![pic](results/2L/GM/3R_Pvalues.png)


#### LEA




#### BAYPASS
