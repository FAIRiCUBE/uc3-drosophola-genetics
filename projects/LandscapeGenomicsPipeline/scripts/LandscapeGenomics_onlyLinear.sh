#!/bin/sh

################################################# CREATING REQUIRED DIRECTORIES TO PROCESS AND STORE OUTPUTS ##################################################
#scriptdir="/home/sonjastndl/s3/LGA/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts"
LGAdir=$(pwd)
scriptdir="${LGAdir}/scripts"
wd="$1"
continent="$2"
arm="$3"
resultsdir="${wd}/results"

mkdir $wd
mkdir $resultsdir
pdir="${resultsdir}/${arm}"
mkdir $pdir

FinalOut="${resultsdir}/${arm}/summary"
mkdir $FinalOut

cd $wd
mkdir data

################################################## DOWNLOAD VCF FROM DEST.bio AND EXTRACT SAMPLENAMES  ##################################################
echo "DOWNLOADING VCF FROM DEST.bio AND EXTRACTING SAMPLENAMES"
#
#wget --tries=inf "http://berglandlab.uvadcos.io/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
#wget "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv"
#wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz


### Extract the information on the available samples
#awk '{FS=","}{if (NR!=1) {print $1}}' dest_v2.samps_3May2024.csv > data/samplenames.csv
#cp dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz data/PoolSeq2024.vcf.gz
#cp  dmel-all-r6.57.gff.gz > data/dmel-all-r6.57.gff.gz
#cp /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/Data/Drosophila/dmel-all-r6.57.gff.gz > data/dmel-all-r6.57.gff.gz
#

module load Tools/vcftools_0.1.13

#################################################################### NAMING THE ANALYSIS ###################################################################################

#input=${wd}/data/PoolSeq2024.vcf.gz
input="$5"
echo $input
##
##
#awk -F, 'NR > 1 && $6 == "Europe" {print $1}' dest_v2.samps_3May2024.csv > ${wd}/data/EuropeSamples.csv
#awk -F, -v var="$continent" 'NR > 1 && $6 == var && $46 == "Pass" {print $1}' dest_v2.samps_3May2024.csv > ${wd}/data/EuropeSamples_Pass.csv
#awk -F "," '$(NF-7) !="Pass" || $(NF-9)<15 {print $1"\t"$(NF-7)"\t"$(NF-9)}' dest_v2.samps_3May2024.csv > ${wd}/data/REMOVE.ids
##
samplelist="$4"
#samplelist="${wd}/data/EuropeSamples_Pass.csv"
#samplelist="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData/data/EuropeSamples_Pass.csv"
#WorldClim="${wd}/data/WorldclimData.csv"
#metadata="${wd}/dest_v2.samps_3May2024.csv"
metadata="$6"
#metadata_new="${wd}/data/metadata_new.csv"
envdata="$7"
#
#
############################### MERGE BIOCLIM VARIABLES (e.g. bio1) WITH OTHER ENV DATA (e.g. NSAT) ###############################################################################
#echo "Merging environmental data sources"
#python ${scriptdir}/MergeData.py --biovariable bio1 bio2 --output $metadata_new --samplenames $samplelist --metadata $metadata --worldclim $WorldClim --climate_extra /home/sonjastndl/s3/ClimateData/TESTMARIA_5.csv --climate_var NSAT
#python3 ${scriptdir}/MergeData.py --biovariable bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 --output $metadata_new --samplenames $samplelist --metadata $metadata --worldclim $WorldClim 
#
#################################### Remove polyploidies,subsample population samples and exlcude all sites with missing data #################################################
#
#
###bash /home/sonjastndl/s3/vcftools.sh $input $samplelist $arm
##
#
#echo "Starting VCFtools" to filter random SNPs
echo $input 
echo $samplelist
echo $arm
#echo "START FILTERING"
#
#echo "WD"
#echo $wd
#
##vcftools --gzvcf $input --keep $samplelist  --minDP 15 --stdout --recode-INFO-all --recode | grep -v "\./\." > ${wd}/results/$arm/Subsampled_$arm.recode2_DP15.vcf.gz
#
#pigz -dc $input | awk '$0~/^#/ || length($5)==1' | vcftools --vcf - \
#        --keep $samplelist \
#        --remove ${wd}/data/REMOVE.ids \
#        --stdout --recode-INFO-all \
#        --recode | grep -v "\./\." | gzip  >  ${wd}/results/$arm/Subsampled_$arm.recode2_DP15.vcf.gz
#
#echo "DONE WITH FILTERING"

#gzip ${wd}/results/$arm/Subsampled_$arm.recode2_DP15.vcf.gz

##
Sub2="${wd}/results/${arm}/Subsampled_${arm}.recode2_DP15.vcf.gz"
Sub3="${wd}/results/${arm}/Subsampled_${arm}.recode3_DP15.vcf.gz"
Sub4="${wd}/results/${arm}/Subsampled_${arm}.final_DP15.vcf.gz"
##
##
#################################################### RANDOMLY PICK n LINES FROM VCF ###########################################################
##
##

#echo "SUBSAMPLING WITH PYTHON"
#echo $Sub2
#
#
#
#python3 ${scriptdir}/SubsampleVCF.py \
#    --input ${Sub2} \
#    --snps all \
#    --output ${Sub3}
#
#module load Tools/bcftools-1.16 
###
####################################################    RUN BCFTOOLS    #######################################################################
#bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $Sub3 | gzip > $Sub4



#echo "GETTING INTRONIC SNPs"
#get intronic SNPS
#gff_file="${wd}/data/dmel-all-r6.57.gff.gz"

#neutralSNPs="${wd}/results/${arm}/Subsampled_${arm}_NeutralSNPS_80.tsv"

#echo "FILTERING INTRONIC SNPs"
#
#python3 ${scriptdir}/IntronicSNPS.py \
# --gff $gff_file \
# --vcf ${Sub4} \
# --target-length 80 \
# --output $neutralSNPs
#
#
#vcftools --gzvcf $Sub4 \
#   --positions $neutralSNPs \
#   --recode --stdout | gzip > ${wd}/results/${arm}/Subsampled_neutral.vcf.gz
#
#################################################### CONVERT TO ALLELE FREQUENCIES ############################################################

#echo "CONVERTING VCF FILE TO ALLELE FREQUENCY FILE"
#
#python3 ${scriptdir}/VCF2AF.py --input ${Sub4} > ${wd}/results/${arm}/Subsampled_${arm}.final_DP15.af
#python3 ${scriptdir}/VCF2AF.py --input ${wd}/results/${arm}/Subsampled_neutral.vcf.gz > ${wd}/results/${arm}/Neutral.final.af
#
#
####annotate SNPs#
#
#echo "ANNOTATING SNPs"
#snpEff="/opt/bioinformatics/snpEff/snpEff.jar"
#java -jar /media/inter/ssteindl/DEST/DESTv1_Pipeline/shell/snpEff/snpEff.jar
#
##/usr/lib/jvm/java-11-openjdk-11.0.22.0.7-2.el8.x86_64/bin/java -jar /media/inter/ssteindl/DEST/DESTv1_Pipeline/shell/snpEff/snpEff.jar
#
#module load Tools/snpEff        
#
#annotated="${wd}/results/${arm}/Subsampled_${arm}.final_DP15.ann.vcf.gz"
#/usr/lib/jvm/java-11-openjdk-11.0.22.0.7-2.el8.x86_64/bin/java -jar $snpEff ann BDGP6.28.99 $Sub4 | gzip >> $annotated
#
#more $annotated | gunzip | awk ' !/^#/ {split($8,a,"|"); print $1 " " $2 " " a[4]}' > ${wd}/results/annotations.txt
#awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' ${wd}/results/${arm}/Subsampled_${arm}.final_DP15.af ${wd}/results/annotations.txt > ${wd}/results/annotated_used_frqs.txt
###
#AF="${wd}/results/${arm}/Subsampled_${arm}.final_DP15.af"
AF="$8"
#
###conda deactivate 

################################################## PERFORM LINEAR REGRESSION #################################################################
echo "PERFORMING LINEAR REGRESSION"

#unmute to run 
Rscript ${scriptdir}/Plot_pvalues.R ${wd} $AF $envdata $arm ${LGAdir}/${FinalOut}


##
##bash ${scriptdir}/transpose_and_split.sh $metadata_new 
#
################################################## PERFORM LFMM 2 (LATENT FACTOR MIXED MODEL) ###################################################
#echo "PERFORMING LFMM"
#
#LeaOut="${wd}/results/${arm}/LEA"
#mkdir $LeaOut
##
#variables=$(head -n 1 "$envdata" | sed 's/\r$//')
##echo $variables
###  #Use a for loop to iterate over words (assuming space-separated words)
#IFS=',' read -ra header_elements <<< "$variables"
#
##Rscript /home/sonjastndl/s3/InstallLea.R
###test
##Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $metadata_new "bio1" $rep
#
##for ((i = 1; i < 3; i++)); do
#for ((i = 1; i < ${#header_elements[@]}; i++)); do
#    element=${header_elements[i]} 
#    echo "$element"
#    #echo $LeaOut
#    #echo $AF
#    rep=1
#    echo $rep
#    # Add your processing here
#    Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $envdata $element $rep
#    #Rscript ${scriptdir}/LEA_RunLFMM2.R $LeaOut $AF $metadata_new "bio1" 1
#    #If needed average the Repetitions
#    #Rscript Rscript ${scriptdir}/LEA_ZPcalc.R $LeaOut $nK $nR $AF $var
#done
##Rscript ${scriptdir}/PlotLEAPValues.r $wd $AF $metadata_new $arm $FinalOut
#
####Rscript ${scriptdir}/ComparePValues.R $AF ${wd}/results/${arm}/GM $LeaOut $FinalOut
#
#### RDA
## Get Intronic SNPs
#
##wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.57.gff.gz
#
##cd .. 
#
##wd="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData"
#####RDA
##AF_file="${wd}/results/fullgenome2/Subsampled_fullgenome2.final_DP15.af"
##metadata="${wd}/dest_v2.samps_3May2024.csv"
#neutral_af="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Neutral.final.af"
#RDA_out="${wd}/results/${arm}/RDA"
#mkdir $RDA_out
##RDA_out="${wd}/results/fullgenome2/RDA_annotations"
##wc_folder="" ##needs to be included in script
##env_data="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/wc2.5"
###
#
#
##Rscript /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/UC3_Landscape_RDA.R $AF $metadata $neutral_af $RDA_out $envdata "Europe" "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/annotationdata/annotations.txt"
#Rscript /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullDataRun/results/AllSNPs/RDA/NewApproachRDA.r $AF $metadata $neutral_af $RDA_out $envdata "Europe" "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/annotationdata/annotations.txt"


###Rscript /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/RDA.R $AF_file $metadata $neutral_af $RDA_out $nc_folder


#for value in {1..20}
#do
#    # Assign the current value to the error variable
#    error=$value
#    echo $error
#    # Run the R script with the given parameters
#    Rscript /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/RDA_MAFtest.R \
#    /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k_newSNPs/SubsampledEurope50k.af \
#    /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/dest_v2.samps_3May2024.csv \
#    /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/matchingEurope_neutral.af \
#    /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA \
#    "Europe" \
#    $error
#done