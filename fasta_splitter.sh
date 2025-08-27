#!/bin/bash

## script for splitting a FASTA file into segments to prevent gene stacking issues with iVar
## requires an input FASTA file and an output path for the cleaned FASTA

VARIANT=$1
ORIGINAL_FASTA=$2

## set output directory
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/segmented"

## create the output directory if it does not exist yet
if [ ! -d ${OUTPUT_DIR} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## create an empty cleaned fasta file
CLEANED_FASTA="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/${VARIANT}_reference_cleaned.fasta"

## clean fasta headers if they are not cleaned already
awk '{if($0 ~ /^>/){split($0,a," "); print a[1]} else {print $0}}' $ORIGINAL_FASTA > $CLEANED_FASTA


## create a separate fasta for each segment
awk -v outdir="$OUTPUT_DIR" '
  /^>/ {
    if (out) close(out)
    seg = $1
    sub(/^>/,"",seg)
    out = outdir "/" seg ".fasta"
  }
  { print > out }
' "$CLEANED_FASTA"
