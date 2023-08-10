## GapFilling project

### Link to ImPoolation Repository for heuristic gap filling

https://github.com/capoony/ImPoolation

### Metadata

Metadata file: `projects/gap_filling/MartinData/data/METADATA_dest_v2.samps_25Feb2023.csv`

### generate full AF dataset

In the following, the script `VCF2AF.py` converts the VCF to an tab-separated matrix file, where rows are genomic positions and columns are populations. The first three columns describe the Chromosome, Position and the two most-common alleles, respectively. For each population and position, this file contains counts for the two most-common alleles separated by a comma. Missing data is indicated by `.,.`.

```
#### 1) convert VCF to AF file format, now retaining the information on read depths 

cd /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/

jobs=100
input=/media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz

gunzip -c $input | parallel \
    -k \
    --pipe \
    -j $jobs \
    --no-notice \
    --cat python3 scripts/VCF2AF.py \
        --input {} \
        | gzip > data/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.af.gz
```

### subset North American and European datasets

In the next step, I isolate all populations from North America and from Europe, respectively. In addition, this step only retains position, that are (still) polymorphic in at least one of the populations in each of the continents.

```bash

cd /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/

for continent in North_America Europe; do

    python scripts/SubsetDataByMeta.py \
        --input data/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.af.gz \
        --meta data/dest_v2.samps_8Jun2023.csv \
        --criteria continent=${continent} \
        --out data/${continent}_PoolSNP.001.50.8Jun2023.norep

done
```

### Gap profiles for all populations

Here, I am counting the proportion of missing positions in the dataset for each population. Importantly, I am distinguishing between single missing sites or consecutive sites (of length *n*) along the genome with missing data. Thus, the output files contain rows for missing sites (gaps) of length 1 to length *n*. Every column shows the proportion relative to the total SNP dataset.

```bash

cd /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/

### first create GapProfile based on "real" positions

for continent in North_America Europe; do

    python scripts/GapProfile.py \
        --input data/${continent}_PoolSNP.001.50.8Jun2023.norep.af.gz \
        --meta data/${continent}_PoolSNP.001.50.8Jun2023.norep.meta \
        > data/${continent}_PoolSNP.001.50.8Jun2023.norep.gapprofile.txt &

done

### then create GapProfile based on row indices

for continent in North_America Europe; do

    python scripts/GapProfile_byrows.py \
        --input data/${continent}_PoolSNP.001.50.8Jun2023.norep.af.gz \
        --meta data/${continent}_PoolSNP.001.50.8Jun2023.norep.meta \
        > data/${continent}_PoolSNP.001.50.8Jun2023.norep.gapprofile_byrows.txt &

done


```
