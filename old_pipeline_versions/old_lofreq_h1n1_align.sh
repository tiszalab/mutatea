#!/bin/bash
set -euo pipefail

## a script for pulling all iav_serotype assigned H1N1 reads (.fastq)
## aligning them to an H1N1 reference with minimap2
## samtools sort and index the bam files
## use loseq to get vcf files
## use varmint to get tsv files

POOLID=$1
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/${POOLID}"

if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference.fasta"
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*H1N1.R1.fastq" )
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/GCA_039113225.1_ASM3911322v1_genomic.gff"


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

		## create directory for filtered vcf files
		mkdir -p ${OUTPUT_DIR}/filtered
		
		## lofreq to see variants from reference genome
		lofreq call -f $REF --no-default-filter --min-cov 1 --max-depth 1000000 --force-overwrite --verbose -o ${OUTPUT_DIR}/${SAMPLE}.raw.vcf $BAM
		
		# created a temporary folder so the filter step can overwrite the vcf
		TMP=$(mktemp -u --suffix=.vcf)
		lofreq filter -i ${OUTPUT_DIR}/${SAMPLE}.raw.vcf --cov-min 1 --af-min 1 --verbose -o "$TMP"
		mv "$TMP" ${OUTPUT_DIR}/filtered/${SAMPLE}.vcf

		# Create a VCF file for the sample
        ##bcftools mpileup -Ou -f $REF -d 1000000 -q 0 -Q 0 -a DP,AD $BAM | bcftools view -Oz -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.vcf.gz
		# Index the VCF file
		##bcftools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.vcf.gz
        
        ## varmint to convert vcf to tsv


		## keep only the non-empty tsv files
		## need to consider that every tsv file has a header line
		##if [ $(wc -l < "${OUTPUT_DIR}/${SAMPLE}.tsv") -le 1 ]; then
		##	rm -f "${OUTPUT_DIR}/${SAMPLE}.tsv"
		##else
	    ##	echo "Variants found for ${SAMPLE}"
		##fi
done
fi