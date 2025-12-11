#!/bin/bash

#SBATCH -p cmmr
#SBATCH -w cmp16

VARIANT=$1   

for d in /gpfs1/projects/Pools/EsViritu/p*/; do
    POOLID=$(basename "$d")
    sbatch \
           /gpfs1/projects/Tisza_Lab/crm_flu_mutatome/flu_mutatome_pipelines/09.09_my_mergedbam.sh "$VARIANT" "$POOLID"
done
