## first make SNPs based on chromatin

## working directory
WD=/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling

## Define windows for gap-filling: 1) Based on chromatin (Euchromatin, Centromer, Telomer); 2) Based on Chromosomal arm (X,2L,2R,3L,3R,Y,4); 3) Based on genomic windows including 5mio positions.
## Each of the output files contains three columns: Type (The WindowID), Chrom (the Chromosome name ), Pos (The position on the corresponding chromosome)
python ${WD}/makeWindows/scripts/splitChromatin.py \
    --input /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/Subsample1.vcf.gz \
    --output ${WD}/makeWindows/data/Europe \
    --chromatin ${WD}/makeWindows/data/chromatin.txt \
    --window 5000000

## Convert the SNPs in the VCF file format to allele frequencies (*_freq.csv.gz) for each sample and genomic position and to weights (*_weight.csv.gz) which corresponds to the sequencing depth at every position
python ${WD}/makeWindows/scripts/vcf2af.py \
    --input /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/Subsample1.vcf.gz \
    --output ${WD}/makeWindows/data/Europe_genomewide
