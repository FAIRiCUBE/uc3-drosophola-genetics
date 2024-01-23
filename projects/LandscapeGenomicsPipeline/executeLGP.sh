## Please note, this is not the latest DEST data set
#cd data
#wget --tries=inf "http://berglandlab.uvadcos.io/gds/dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gdsdest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz"
#wget "https://github.com/DEST-bio/DESTv2/blob/main/populationInfo/dest_v2.samps_25Feb2023.csv"
###take all samples
#awk '{FS=","}{if (NR!=1) {print $1}}' $samplefile > samplenames.csv # make this a path
#mv data/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz data/PoolSeq2023.vcf.gz

###extract variable from worldclim file 


bash /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/run_LGP_jobs.sh \
  -d /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline \
  -o NSATvsBIO1_7 \
  -s /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/european.csv\
  -i /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/PoolSeq2023.vcf.gz \
  -m /media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest_v2.samps_26April2023.csv \
  -w /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/dest.worldclim.csv \
  -b bio1 \
  -v /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv 


bash /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/run_LGP_jobs.sh \
  -d /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline \
  -o NSATvsBIO1_8 \
  -s /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/cline8.csv\
  -i /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/PoolSeq2023.vcf.gz \
  -m /media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest_v2.samps_26April2023.csv \
  -w /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/dest.worldclim.csv \
  -b bio1 \
  -v /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv 


bash /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/scripts/runLGPgenomewide.sh \
  -d /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline \
  -o NSATvsBIO1_9 \
  -s /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/cline8.csv\
  -i /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/PoolSeq2023.vcf.gz \
  -m /media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest_v2.samps_26April2023.csv \
  -w /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/dest.worldclim.csv \
  -b bio1 \
  -v /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv 
