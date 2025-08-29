#!/bin/bash

## this script is for merging BAM files by month and region based on metadata file
## needs to be run after 08.27_varmintonly.sh script has been run to create the BAM files

## set variable
VARIANT=$1

## set output directory and path to metadata file
BASE_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align"
OUTPUT_DIR="${BASE_DIR}/bam_merger_output"
OUTPUT_BAM_DIR="${OUTPUT_DIR}/merged_bams"
META_FILE="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/metadata_combined.csv"

## create output directories if they do not already exist
if [ ! -d "${OUTPUT_DIR}" ] ; then
	mkdir -p "${OUTPUT_DIR}"
fi
if [ ! -d "${OUTPUT_BAM_DIR}" ] ; then
	mkdir -p "${OUTPUT_BAM_DIR}"
fi

# column numbers from metadata CSV
SAMPLE_COL=4   ## column for sample name
DATE_COL=11    ## column for Month.year in the format of 08.2024
REGION_COL=10  ## column for region
POOLID_COL=8   ## column for poolID

## create a list file for each month-region combination
while IFS=',' read -r -a cols; do
    SAMPLE=$(echo "${cols[$SAMPLE_COL]}" | xargs)  # Trim whitespace
    DATE=$(echo "${cols[$DATE_COL]}" | xargs)
    REGION=$(echo "${cols[$REGION_COL]}" | xargs)
    POOLID=$(echo "${cols[$POOLID_COL]}" | xargs)

    ## set path to BAM files
    BAM_PATH="${BASE_DIR}/${POOLID}/${SAMPLE}.${POOLID}.sort.bam"

    ## want each list file to be named by the month and region it is being grouped by
    MERGE_KEY="${DATE}.${REGION}"
    LIST_FILE="${OUTPUT_DIR}/${MERGE_KEY}.list"

    if [[ -f "${BAM_PATH}" ]]; then
        echo "${BAM_PATH}" >> "${LIST_FILE}"
    fi
done < <(tail -n +2 "${META_FILE}")  # skip header

## merge BAMs for each month-region combination
for LIST_FILE in "${OUTPUT_DIR}"/*.list; do
    if [[ -f "${LIST_FILE}" && -s "${LIST_FILE}" ]]; then  # ensure the list file exists and is not empty
        MERGE_KEY=$(basename "${LIST_FILE}" .list)
        MERGED_BAM="${OUTPUT_BAM_DIR}/${MERGE_KEY}.bam"

        echo "Merging the BAMs for ${MERGE_KEY}"
        samtools merge -f "${MERGED_BAM}" -b "${LIST_FILE}"
    else
        echo "WARNING: Skipping empty or missing list file: ${LIST_FILE}" >&2
    fi
done