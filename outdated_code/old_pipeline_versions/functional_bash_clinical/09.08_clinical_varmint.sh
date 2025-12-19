#!/bin/bash

## script for using varmint on sorted and indexed merged BAM outputs to capture variants in tsv format
## run after running 09.08_clinical.sh to create, sort, and index the merged BAM files

## set variable
VARIANT=$1

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## set directories
BASE_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/${VARIANT}"
GROUPED_BAMS="${BASE_DIR}/grouped_bams/"
TSV_OUTPUT="${BASE_DIR}/tsv_output"
META_FILE="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/metadata/${VARIANT}_clinical_md_my.tsv"

## create the output directory if it does not already exist
if [ ! -d ${TSV_OUTPUT} ] ; then 
    mkdir -p ${TSV_OUTPUT}
fi

## load in reference files
REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/${VARIANT}_reference.fasta"
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/gff/${VARIANT}.gff"

## run varmint on the grouped BAM files
for GROUPED_BAM in "${GROUPED_BAMS}"/*.sort.bam; do
    if [[ -f "${GROUPED_BAM}" ]]; then  ## ensure the sorted BAM file exists
        SAMPLE=$(basename "${GROUPED_BAM}" .sort.bam)
        TSV_OUTPUT_FILE="${TSV_OUTPUT}/${VARIANT}.${SAMPLE}.tsv"

        ## run varmint
        varmint --bam "${GROUPED_BAM}" --ref "${REF}" --gff "${GFF}" --out "${TSV_OUTPUT_FILE}"

        ## check if TSV file exists and is not empty
        if [[ -f "${TSV_OUTPUT_FILE}" && $(wc -l < "${TSV_OUTPUT_FILE}") -le 1 ]]; then
            rm -f "${TSV_OUTPUT_FILE}"
        elif [[ -f "${TSV_OUTPUT_FILE}" ]]; then
            echo "Variants found for ${SAMPLE}"
        fi
    fi
done