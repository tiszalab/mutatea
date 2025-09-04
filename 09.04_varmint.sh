#!/bin/bash

## script for using varmint on merged BAM outputs to capture variants in tsv format

## set variables
VARIANT=$1
POOLID=$2

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## set directories
BASE_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align"
OUTPUT_DIR="${BASE_DIR}/pools/${POOLID}"
MERGED_BAM_DIR="${BASE_DIR}/bam_merger_output/merged_bams"
TSV_OUTPUT="${BASE_DIR}/bam_merger_output/tsv_files" 
META_FILE="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/metadata_combined.csv"

## create the output directories if they do not already exist
if [ ! -d ${OUTPUT_DIR} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi
if [ ! -d ${MERGED_BAM_DIR} ] ; then
	mkdir -p ${MERGED_BAM_DIR}
fi
if [ ! -d ${TSV_OUTPUT} ] ; then 
    mkdir -p ${TSV_OUTPUT}
fi

## run varmint on the merged BAM files
for SORTED_BAM in "${MERGED_BAM_DIR}"/*.sort.bam; do
    if [[ -f "${SORTED_BAM}" ]]; then  ## ensure the sorted BAM file exists
                BAM="${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"
        SAMPLE=$(basename "${SORTED_BAM}" .sort.bam)
        TSV_OUTPUT_FILE="${TSV_OUTPUT}/${SAMPLE}.tsv"

        ## run varmint
        varmint --bam "${SORTED_BAM}" --ref /gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/FASTA/H1N1_reference_cleaned.fasta --gff /gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/H1N1.gff --out "${TSV_OUTPUT_FILE}"

        ## check if TSV file exists and is not empty
        if [[ -f "${TSV_OUTPUT_FILE}" && $(wc -l < "${TSV_OUTPUT_FILE}") -le 1 ]]; then
            rm -f "${TSV_OUTPUT_FILE}"
        elif [[ -f "${TSV_OUTPUT_FILE}" ]]; then
            echo "Variants found for ${SAMPLE}"
        fi
    fi
done