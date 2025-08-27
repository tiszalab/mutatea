#!/bin/bash

## script for splitting a GFF file into non-overlapping segments to prevent gene stacking issues with iVar
## requires an input GFF file and an output path for the cleaned GFF

OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/segmented"

if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in file
ORIGINAL_GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/H1N1.gff"

## keep headers from original GFF
HEADERS=$(grep "^#" "$ORIGINAL_GFF")

## get list of all segments
SEGMENTS=$(grep -v "^#" "$ORIGINAL_GFF" | cut -f1 | sort -u)

## create a GFF file for each segment
for SEG in $SEGMENTS; do
    SEGMENTED_GFF="${OUTPUT_DIR}/${SEG}.gff"
    
    # write headers
    echo "$HEADERS" > "$SEGMENTED_GFF"
    
    # write only lines belonging to this segment
    awk -v seg="$SEG" '$1==seg {print}' "$ORIGINAL_GFF" >> "$SEGMENTED_GFF"
    
    echo "Created $SEGMENTED_GFF"
done
