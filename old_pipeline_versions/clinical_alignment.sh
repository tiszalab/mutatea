#!/bin/bash

## this script creates BAM files for each sample in the split FASTA files
## this script should be run after accession_fasta_splitter.sh to align the split FASTA files to the reference genome

## set variable
VARIANT=$1

## set output directory
INPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/${VARIANT}_split_fasta_files"
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/${VARIANT}_bam_files"

## create output directory if it does not already exist
if [ ! -d "${OUTPUT_DIR}" ] ; then
	mkdir -p "${OUTPUT_DIR}"
fi

## set input and output files
REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/${VARIANT}_reference_cleaned.fasta"

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome


for FASTA in "${INPUT_DIR}"/*.fasta; do
    ACC=$(basename "$FASTA" .fasta)

    ## minimap2 alignment
    minimap2 -t 48 -a "${REF}" "${FASTA}" > "${OUTPUT_DIR}/${ACC}.sam"

    ## samtools convert to BAM, sort, and index
    samtools view -@ 48 -bS "${OUTPUT_DIR}/${ACC}.sam" > "${OUTPUT_DIR}/${ACC}.bam"
    samtools sort -@ 48 -o "${OUTPUT_DIR}/${ACC}.sort.bam" "${OUTPUT_DIR}/${ACC}.bam"
    samtools index "${OUTPUT_DIR}/${ACC}.sort.bam"
done

echo "Alignments done for all ${VARIANT} clinical samples"