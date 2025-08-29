#!/bin/bash

## a script for pulling all iav_serotype assigned H1N1 reads (.fastq)
## aligning them to an H1N1 reference with minimap2
## samtools sort and index the bam files
## use ivar to get tsv files
## should run the fasta_splitter.sh script and gff_splitter.sh first to create segmented FASTA and GFF files

## set variables
VARIANT=$1
POOLID=$2

## set output directories
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/${POOLID}"
TSV_OUTPUT="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/${POOLID}/tsv_files"

## create the output directories if they do not already exist
if [ ! -d ${OUTPUT_DIR} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi
if [ ! -d ${TSV_OUTPUT} ] ; then
    mkdir -p ${TSV_OUTPUT}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in reference
REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/FASTA/${VARIANT}_reference_cleaned.fasta"
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/gff/${VARIANT}.gff"

## load in reads
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*${VARIANT}.R1.fastq" )

if [ ! -z "$R1_LIST" ] ; then

	echo "$R1_LIST" | while read READ1 ; do

		echo $READ1
		date

		READ2=$( echo $READ1 | sed 's/R1.fastq/R2.fastq/g' )

		R1_BASE=$( basename $READ1 )

		SAMPLE=$( echo $R1_BASE | cut -d "." -f 1 )

		## minimap2
        minimap2 -ax sr $REF $READ1 $READ2 > ${OUTPUT_DIR}/${SAMPLE}.sam 

        ## samtools sort and index
        samtools view -@ 48 -bS  ${OUTPUT_DIR}/${SAMPLE}.sam > ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
		samtools sort -@ 48 -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
        samtools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam

        BAM="${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"

		varmint --bam $BAM --ref $REF --gff $GFF --out "${TSV_OUTPUT}/${SAMPLE}.tsv"

        if [ $(wc -l < "${TSV_OUTPUT}/${SAMPLE}.tsv") -le 1 ]; then
            rm -f "${TSV_OUTPUT}/${SAMPLE}.tsv"
        else
            echo "Variants found for ${SAMPLE}"
        fi
    done
fi