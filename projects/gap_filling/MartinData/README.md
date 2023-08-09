## GapFilling project

### Link to ImPoolation Repository for heuristic gap filling

https://github.com/capoony/ImPoolation

### Metadata

Metadata file: `projects/gap_filling/MartinData/data/METADATA_dest_v2.samps_25Feb2023.csv`

### generate full AF dataset

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

```bash

cd /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/

for continent in North_America Europe; do

    python scripts/GapProfile.py \
        --input data/${continent}_PoolSNP.001.50.8Jun2023.norep.af.gz \
        --meta data/${continent}_PoolSNP.001.50.8Jun2023.norep.meta \
        > data/${continent}_PoolSNP.001.50.8Jun2023.norep.gapprofile.txt

done
```
