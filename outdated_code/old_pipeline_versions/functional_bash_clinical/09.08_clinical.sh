#!/bin/bash

## this script extracts clinical flu samples from a merged FASTA, grouping them based on the Month_Year column from their respective metadata

## set variables
VARIANT=$1

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in files
METADATA="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/metadata/${VARIANT}_clinical_md_my.tsv"
FASTA="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/${VARIANT}_clinical.fasta"
REF_FASTA="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/${VARIANT}_reference.fasta"

## set output directory
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/clinical_flu_mutatome/${VARIANT}/grouped_bams"

## create output directory if it does not already exist
if [ ! -d ${OUTPUT_DIR} ] ; then 
    mkdir -p ${OUTPUT_DIR}
fi

## get unique values from the Month_Year column of metadata
month_years=$(tail -n +2 "$METADATA" | cut -f15 | sort | uniq)


## loop for each unique Month_Year value
for month_year in $month_years; do
    echo "Processing $month_year..."

    ## extract sequence IDs for the current month_year group
    seq_ids=$(awk -v my="$month_year" -F'\t' '$15 == my {print $1}' "$METADATA")

    ## create a FASTA file for each month_year group
    if [ -n "$seq_ids" ]; then
        ## Log the sequence IDs being merged
        SEQ_LOG="${OUTPUT_DIR}/${month_year}_seq_ids.txt"
        echo "$seq_ids" > "$SEQ_LOG"

        echo "$seq_ids" | awk '
            NR==FNR { ids[$1]; next }   # first pass: read IDs into array
            /^>/ {
                header=$0
                split($1,a,"|")         # take first token after ">"
                split(a[1],b," ")
                id=substr(b[1],2)       # remove ">"
            }
            { if (id in ids) print }    # print header+sequence if ID matches
        ' - "$FASTA" > "${OUTPUT_DIR}/${month_year}.fasta"

        ## check if the FASTA file is still empty
        if [ ! -s "${OUTPUT_DIR}/${month_year}.fasta" ]; then
            echo "Warning: FASTA file for $month_year is empty. Check metadata and input FASTA."
        else
            ## align reads with minimap2 with asm5 preset
            minimap2 -ax asm5 "$REF_FASTA" "${OUTPUT_DIR}/${month_year}.fasta" > "${OUTPUT_DIR}/${month_year}.sam"

            ## sort and index BAM
            samtools sort -o "${OUTPUT_DIR}/${month_year}.sort.bam" "${OUTPUT_DIR}/${month_year}.sam"
            samtools index "${OUTPUT_DIR}/${month_year}.sort.bam"
        fi
    else
        echo "No ${VARIANT} sequences found for $month_year"
    fi
done
