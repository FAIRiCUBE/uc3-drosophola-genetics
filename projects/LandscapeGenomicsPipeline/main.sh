##!/bin/sh
## name of Job
#PBS -N RasdaLEA_only
## Redirect output stream to this file.
#PBS -o /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RasdaTest/logs
## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe
## Select a maximum of 20 cores and 200gb of RAM
#PBS -l select=1:ncpus=8:mem=50gb
## load all necessary software into environment


scriptdir="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts"

cd $scriptdir 

wd="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullDataRun"
continent="Europe" 
arm="AllSNPs"
samplelist="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/data/EuropeSamples_Pass.csv"
input="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/data/PoolSeq2024.vcf.gz"
metadata="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/dest_v2.samps_3May2024.csv"
envdata="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/use_case/rasdaman_data.csv"
AF="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome.final_DP15.af"


#bash LandscapeGenomics.sh $wd $continent $arm $samplelist $input $metadata $envdata $AF >> FullDataRun_LEA_qsub.log 2>&1
#bash LandscapeGenomics_qsub.sh $wd $continent $arm $samplelist $input $metadata $envdata $AF >> FullDataRun_LEA_qsub.log 2>&1
bash ${scriptdir}/LandscapeGenomics_onlyRDA.sh $wd $continent $arm $samplelist $input $metadata $envdata $AF >> FullDataRun_RDA_qsub.log 2>&1
