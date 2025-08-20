#!/bin/bash

## a script for pulling all iav_serotype assigned H1N1 reads (.fastq)
## aligning them to an H1N1 reference with minimap2
## samtools sort and index the bam files
## use bcftools to get vcf files

POOLID=$1
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/${POOLID}"


if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference_clean.fasta"
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*H1N1.R1.fastq" )

if [ ! -z "$R1_LIST" ] ; then

	echo "$R1_LIST" | while read READ1 ; do

		echo $READ1
		date

		READ2=$( echo $READ1 | sed 's/R1.fastq/R2.fastq/g' )

		R1_BASE=$( basename $READ1 )

		SAMPLE=$( echo $R1_BASE | cut -d "." -f 1 )

		# minimap2
        minimap2 -ax sr $REF $READ1 $READ2 > ${OUTPUT_DIR}/${SAMPLE}.sam 

        # samtools sort and index
        samtools view -@ 48 -bS  ${OUTPUT_DIR}/${SAMPLE}.sam > ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
		samtools sort -@ 48 -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
        samtools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam

        # Set bam variable
        BAM="${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"

        # Create a VCF file for the sample
        bcftools mpileup -Ou -f $REF -q 0 -Q 0 --min-BQ 0 --max-depth 1000000 -a DP,AD $BAM | bcftools view -Oz -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.vcf.gz
        
        #bcftools norm -Oz -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.vcf
        
        # crm set to DNR Index the VCF file
        # bcftools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.vcf.gz
done
fi