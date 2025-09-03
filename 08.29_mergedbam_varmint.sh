#!/bin/bash

## script for aligning reads, merging BAM files, and varmint to capture variants

## set variables
VARIANT=$1
POOLID=$2

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

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in reference files
REF="${BASE_DIR}/reference/FASTA/${VARIANT}_reference_cleaned.fasta"
GFF="${BASE_DIR}/reference/gff/${VARIANT}.gff"

## create BAM files
R1_LIST=$(find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*${VARIANT}.R1.fastq")

if [ ! -z "$R1_LIST" ]; then
    echo "$R1_LIST" | while read READ1; do
        READ2=$(echo "$READ1" | sed 's/R1.fastq/R2.fastq/g')
        R1_BASE=$(basename "$READ1")
        SAMPLE=$(echo "$R1_BASE" | cut -d "." -f 1)

        ## minimap2
        minimap2 -ax sr "$REF" "$READ1" "$READ2" > "${OUTPUT_DIR}/${SAMPLE}.sam"

        ## samtools sort and index
        samtools view -@ 48 -bS "${OUTPUT_DIR}/${SAMPLE}.sam" > "${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam"
        samtools sort -@ 48 -o "${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam" "${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam"
        samtools index "${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"
    done
else 
    ## remove empty POOLID folders
    echo "No R1 files found for POOLID: ${POOLID}"
    rmdir "${OUTPUT_DIR}" 2>dev/null || true
fi


## merge BAM files by metadata
## get column numbers from metadata CSV
SAMPLE_COL=4   ## column for sample name
DATE_COL=11    ## column for Month.year in the format of 08.2024
REGION_COL=10  ## column for region
POOLID_COL=8   ## column for poolID

while IFS=',' read -r -a cols; do
    SAMPLE=$(echo "${cols[$SAMPLE_COL]}" | xargs)
    DATE=$(echo "${cols[$DATE_COL]}" | xargs)
    REGION=$(echo "${cols[$REGION_COL]}" | xargs)
    POOLID=$(echo "${cols[$POOLID_COL]}" | xargs)

    ## set path to BAM files
    BAM_PATH="${BASE_DIR}/pools/${POOLID}/${SAMPLE}.${POOLID}.sort.bam"

    ## want each list file to be named by the month and region combination it is being grouped by
    MERGE_KEY="${DATE}.${REGION}"
    LIST_FILE="${BASE_DIR}/bam_merger_output/${MERGE_KEY}.list"

    if [[ -f "${BAM_PATH}" ]]; then
        echo "${BAM_PATH}" >> "${LIST_FILE}"
        sort -u "${LIST_FILE}" -o "${LIST_FILE}"  ## remove duplicate entries
    fi
done < <(tail -n +2 "${META_FILE}") ## skip header

## merge BAMs for each month-region combination
for LIST_FILE in "${BASE_DIR}/bam_merger_output"/*.list; do
    if [[ -f "${LIST_FILE}" ]]; then
        MERGE_KEY=$(basename "${LIST_FILE}" .list)
        MERGED_BAM="${MERGED_BAM_DIR}/${MERGE_KEY}.bam"

        echo "Merging BAMs for ${MERGE_KEY}"
        samtools merge -f "${MERGED_BAM}" -b "${LIST_FILE}"

        ## sort and index the merged BAM files
        SORTED_BAM="${MERGED_BAM%.bam}.sort.bam"  ## set the variable for the sorted BAM file path
        samtools sort -@ 48 -o "${SORTED_BAM}" "${MERGED_BAM}"
        samtools index "${SORTED_BAM}"
    fi
done

## run varmint on the merged BAM files
for SORTED_BAM in "${MERGED_BAM_DIR}"/*.sort.bam; do
    if [[ -f "${SORTED_BAM}" ]]; then  ## ensure the sorted BAM file exists
        SAMPLE=$(basename "${SORTED_BAM}" .sort.bam)
        TSV_OUTPUT_FILE="${TSV_OUTPUT}/${SAMPLE}.tsv"

        ## run varmint
        varmint --bam "${SORTED_BAM}" --ref "${REF}" --gff "${GFF}" --out "${TSV_OUTPUT_FILE}"

        ## check if TSV file exists and is not empty
        if [[ -f "${TSV_OUTPUT_FILE}" && $(wc -l < "${TSV_OUTPUT_FILE}") -le 1 ]]; then
            rm -f "${TSV_OUTPUT_FILE}"
        elif [[ -f "${TSV_OUTPUT_FILE}" ]]; then
            echo "Variants found for ${SAMPLE}"
        fi
    fi
done
