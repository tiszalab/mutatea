#!/bin/bash

## a script for pulling all iav_serotype assigned H1N1 reads (.fastq)
## aligning them to an H1N1 reference with minimap2
## samtools sort and index the bam files
## use ivar to get tsv files

POOLID=$1
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/${POOLID}"
TSV_OUTPUT="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/${POOLID}/tsv_files"

if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi
if [ ! -d ${TSV_OUTPUT} ] ; then
    mkdir -p ${TSV_OUTPUT}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference_clean.fasta"

## directories for segmented GFF and FASTA files
REF_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/FASTA/segmented"
GFF_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/segmented"

## load in reads
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*H1N1.R1.fastq" )

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

		## determine which segment this bam uses and use the respective FASTA and GFF files
		 samtools idxstats $BAM | awk '$3 > 0 {print $1}' | while read SEG ; do
            SEG_REF="${REF_DIR}/${SEG}.fasta"
            SEG_GFF="${GFF_DIR}/${SEG}.gff"

            if [[ -f $SEG_REF && -f $SEG_GFF ]]; then
                samtools mpileup -A -d 0 -B -Q 0 -f $SEG_REF $BAM | \
                ivar variants -p ${TSV_OUTPUT}/${SAMPLE} -t 0 -m 1 -r $SEG_REF -g $SEG_GFF

                if [ $(wc -l < "${TSV_OUTPUT}/${SAMPLE}.tsv") -le 1 ]; then
                    rm -f "${TSV_OUTPUT}/${SAMPLE}.tsv"
                else
                    echo "Variants found for ${SAMPLE} on $SEG"
                fi
            fi
        done

    done
fi