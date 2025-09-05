#!/bin/bash

## this script splits the combined FASTA file of complete Influenza A genomes into separate FASTA files based on their accession numbers

## set variable
VARIANT=$1

## set output directory
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome"/${VARIANT}_split_fasta_files

## create output directory if it does not already exist
if [ ! -d "${OUTPUT_DIR}" ] ; then
	mkdir -p "${OUTPUT_DIR}"
fi

## set input and output files
COMBINED_FASTA="${VARIANT}_clinical.fasta"

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## split combined FASTA file into separate FASTA files based on accession numbers
awk -v outdir="$OUTPUT_DIR" '
  /^>/ {
    # extract strain name inside parentheses: e.g. (A/British_Columbia/PHL-2032-recombinant/2025(H5N1))
    match($0, /\(A\/[^)]+\)/, arr)
    if (arr[0] != "") {
      strain = arr[0]
      gsub(/[()]/,"",strain)      # remove parentheses
      gsub(/\//,"_",strain)       # replace slashes with _
    } else {
      strain = substr($1,2)       # fallback to accession
    }
    file = outdir "/" strain ".fasta"
  }
  { print >> file }
' "$COMBINED_FASTA"