#!/bin/bash

VARIANT="H1N1"   

for d in /gpfs1/projects/Pools/EsViritu/p*/; do
    POOLID=$(basename "$d")
    sbatch --partition=gpu \
           /gpfs1/projects/Tisza_Lab/crm_flu_mutatome/flu_mutatome_pipelines/08.29_mergedbam_varmint.sh "$VARIANT" "$POOLID"
done
