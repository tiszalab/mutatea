#!/bin/bash

## script to sort and index merged BAM files

## load directories and variables
MERGED_BAM_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/bam_merger_output/merged_bams"

## sort and index the merged BAM files
for MERGED_BAM in "${MERGED_BAM_DIR}"/*.bam; do
    if [[ -f "${MERGED_BAM}" ]]; then   ## ensure the BAM file exists
        SORTED_BAM="${MERGED_BAM%.bam}.sort.bam"  ## set the variable for the sorted BAM file path
        samtools sort -@ 48 -o "${SORTED_BAM}" "${MERGED_BAM}"
        samtools index "${SORTED_BAM}"
    fi
done

