#!/bin/bash
set -euo pipefail

## Script for processing H1N1 samples:
## - Align reads with minimap2 (short-read preset)
## - Sort and index BAM
## - Rescue soft-clipped reads
## - Call variants with iVar

POOLID=$1
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/${POOLID}"

mkdir -p ${OUTPUT_DIR}

# Activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference.fasta"
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/GCA_039113225.1_ASM3911322v1_genomic.gff"

R1_LIST=$(find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*H1N1.R1.fastq")

if [ ! -z "$R1_LIST" ]; then
    echo "$R1_LIST" | while read READ1; do
        echo "Processing $READ1"
        date

        READ2=$(echo $READ1 | sed 's/R1.fastq/R2.fastq/g')
        SAMPLE=$(basename $READ1 | cut -d "." -f1)

        ## Align reads
        minimap2 -ax sr $REF $READ1 $READ2 > ${OUTPUT_DIR}/${SAMPLE}.sam

        ## Convert SAM -> BAM, sort, index
        samtools view -@ 48 -bS ${OUTPUT_DIR}/${SAMPLE}.sam > ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
        samtools sort -@ 48 -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
        samtools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam

        BAM="${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"

        ## Rescue soft-clipped bases
        samtools view -h $BAM | \
        awk 'BEGIN{OFS="\t"} /^@/ {print; next} {gsub(/[0-9]+S/, "&"); print}' | \
        samtools view -bS - > ${BAM%.bam}.softclip.bam

        samtools sort -@ 48 -o ${BAM%.bam}.softclip.sort.bam ${BAM%.bam}.softclip.bam
        samtools index ${BAM%.bam}.softclip.sort.bam

        BAM_SOFT="${BAM%.bam}.softclip.sort.bam"

        ## Call variants with iVar
        samtools mpileup -A -d 0 -B -Q 0 -f $REF $BAM_SOFT | \
            ivar variants -p ${OUTPUT_DIR}/${SAMPLE}_variants -t 0 -m 1 -r $REF -g $GFF

        echo "Finished $SAMPLE"
    done
fi
