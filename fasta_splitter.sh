#!/bin/bash

## script for splitting a FASTA file into segments to prevent gene stacking issues with iVar
## requires an input FASTA file and an output path for the cleaned FASTA

OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/FASTA/segmented"

if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in file
ORIGINAL_FASTA="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference_clean.fasta"

awk '/^>/{ 
    # close previous file
    if (out) close(out)
    # create new file using first word of header (without ">")
    seg=$1
    sub(/^>/,"",seg)
    out=sprintf("%s/%s.fasta", "'$OUTPUT_DIR'", seg)
}
{ print > out }' "$ORIGINAL_FASTA"