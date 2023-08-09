## GapFilling project

### Link to ImPoolation Repository for heuristic gap filling

https://github.com/capoony/ImPoolation

### Metadata

Metadata file: `projects/gap_filling/MartinData/data/METADATA_dest_v2.samps_25Feb2023.csv`

### generate full AF dataset

```
#### 1) convert VCF to AF file format, now retaining the information on read depths 

jobs=100
input=/media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz

gunzip -c $input | parallel \
    -k \
    --pipe \
    -j $jobs \
    --no-notice \
    --cat python3 /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/scripts/VCF2AF.py \
        --input {} \
        | gzip > /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/data/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.af.gz
```

### subset North American dataset

```bash

python /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/scripts/SubsetDataByMeta.py \
    --input /media/inter/mkapun/projects/uc3-drosophola-genetics/projects/gap_filling/MartinData/data/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.af.gz \
    --meta D:\GitHub\uc3-drosophola-genetics\projects\gap_filling\MartinData\data\dest_v2.samps_8Jun2023.csv \
    --criteria continent=North_America \
    --out
    
```