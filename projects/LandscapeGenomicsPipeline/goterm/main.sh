WD=/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/goterm
DIR=/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/goterm/data/LGA_Linear
cd ${WD}/scripts

# get GOWINDA
#curl -O "https://master.dl.sourceforge.net/project/gowinda/Gowinda-1.12.jar?viasf=1"

cd ${WD}/data

# gene annotations
wget ftp://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz
# GO descriptions
curl -O https://current.geneontology.org/ontology/subsets/goslim_agr.obo
# Gene models
curl -O http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gtf/dmel-all-r6.59.gtf.gz
gunzip ${WD}/data/dmel-all-r6.59.gtf.gz

# make associations
python ${WD}/scripts/scripts/MakeGOassociations.py \
    --GODescriptions ${WD}/data/goslim_agr.obo \
    --GOAnnotations ${WD}/data/gene_association.fb.gz \
    > geneset.tsv

java -jar ${WD}/scripts/scripts/Gowinda-1.12.jar \
    --annotation-file ${WD}/data/dmel-all-r6.59.gtf \
    --gene-set-file ${WD}/data/geneset.tsv \
    --snp-file ${WD}/data/all_loci.tsv \
    --candidate-snp-file /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/PermutationOutliers2.tsv \
    --simulations 1000000 \
    --threads 100 \
    --mode snp \
    --output-file /media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/goterm/results/RDA/AllVarsPermuted


#for file in ${DIR}/*.tsv
#do
#echo $file
#prefix=$(basename "$file" | sed 's/\..*//')
#echo "$prefix"
#java -jar ${WD}/scripts/scripts/Gowinda-1.12.jar \
#    --annotation-file ${WD}/data/dmel-all-r6.59.gtf \
#    --gene-set-file ${WD}/data/geneset.tsv \
#    --snp-file ${WD}/data/all_loci.tsv \
#    --candidate-snp-file ${file} \
#    --simulations 1000000 \
#    --threads 100 \
#    --mode snp \
#    --output-file ${WD}/results/Linear_approach/${prefix}
#done