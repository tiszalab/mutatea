#!/bin/bash

## a script for pulling all iav_serotype assigned H1N1 reads (.fastq)
## aligning them to an H1N1 reference with minimap2
## samtools sort and index the bam files
## use bcftools to get vcf files

POOLID=$1


if [ ! -d ${POOLID} ] ; then
	mkdir ${POOLID}
fi

source /cmmr/prod/envParams/condanewenv.init && conda activate crm_pangenome

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference.fasta"
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*H1N1.R1.fastq" )


if [ ! -z "$R1_LIST" ] ; then

	echo "$R1_LIST" | while read READ1 ; do

		echo $READ1
		date

		READ2=$( echo $READ1 | sed 's/R1.fastq/R2.fastq/g' )

		R1_BASE=$( basename $READ1 )

		SAMPLE=$( echo $R1_BASE | cut -d "." -f 1 )

		## minimap2
        /gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/minimap2-2.30_x64-linux/minimap2 -ax sr $REF $READ1 $READ2 > ${TMPDIR}/${SAMPLE}.sam 

        ## samtools sort and index
        samtools view -@ 48 -bS  ${TMPDIR}/${SAMPLE}.sam > ${TMPDIR}/${SAMPLE}.${POOLID}.bam
		samtools sort -@ 48 -o ${POOLID}/${SAMPLE}.${POOLID}.sort.bam ${TMPDIR}/${SAMPLE}.${POOLID}.bam
        samtools index ${POOLID}/${SAMPLE}.${POOLID}.sort.bam

        ## bcftools to get vcf
        BAM="${POOLID}/${SAMPLE}.${POOLID}.sort.bam"
        # Create a VCF file for the sample
        bcftools mpileup -Ou -f $REF $BAM | bcftools call -mv -Oz -o ${POOLID}/${SAMPLE}.vcf.gz
        
        # Index the VCF file
        bcftools index ${POOLID}/${SAMPLE}.vcf.gz
done
fi