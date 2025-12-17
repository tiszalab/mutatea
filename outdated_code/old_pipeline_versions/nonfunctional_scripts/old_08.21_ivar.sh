#!/bin/bash

## influenza variant calling with a cleaned GFF to prevent gene stacking issue with iVar

VARIANT=$1
POOLID=$2
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/${POOLID}"

if [ ! -d ${POOLID} ] ; then
	mkdir -p ${OUTPUT_DIR}
fi

## activate conda environment
source /cmmr/prod/envParams/condanewenv.init && conda activate crm_flutatome

## load in files
REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/${VARIANT}_reference.fasta"
R1_LIST=$( find /gpfs1/projects/Pools/EsViritu/${POOLID} -type f -name "*${VARIANT}.R1.fastq" )
ORIGINAL_GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/${VARIANT}_align/reference/gff/${VARIANT}.gff"

## create an empty file to catch output
GFF_CLEAN="${OUTPUT_DIR}/non_overlapping.gff"

## split the GFF into non-overlapping segments by splitting overlapping regions
awk 'BEGIN{OFS="\t"}
/^#/ {print; next}
$3 == "CDS" {
    # Extract gene name
    if (match($9, /gene=([^;]+)/, arr)) {
        gene = arr[1]
    } else if (match($9, /Name=([^;]+)/, arr)) {
        gene = arr[1]
    } else {
        gene = "unknown"
    }
    
    # Store all features
    features[NR] = $0
    starts[NR] = $4
    ends[NR] = $5
    genes[NR] = gene
    regions[NR] = $1
    strands[NR] = $7
    attrs[NR] = $9
}
END {
    # Sort by start position
    n = NR
    for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (starts[i] > starts[j]) {
                # Swap
                temp = starts[i]; starts[i] = starts[j]; starts[j] = temp
                temp = ends[i]; ends[i] = ends[j]; ends[j] = temp
                temp = genes[i]; genes[i] = genes[j]; genes[j] = temp
                temp = regions[i]; regions[i] = regions[j]; regions[j] = temp
                temp = strands[i]; strands[i] = strands[j]; strands[j] = temp
                temp = attrs[i]; attrs[i] = attrs[j]; attrs[j] = temp
            }
        }
    }
    
    # Create non-overlapping segments
    prev_end = 0
    for (i = 1; i <= n; i++) {
        if (starts[i] > prev_end) {
            # Non-overlapping, use as is
            print regions[i], ".", "CDS", starts[i], ends[i], ".", strands[i], ".", attrs[i]
            prev_end = ends[i]
        } else {
            # Overlapping, adjust start position
            new_start = prev_end + 1
            if (new_start <= ends[i]) {
                print regions[i], ".", "CDS", new_start, ends[i], ".", strands[i], ".", attrs[i]
                prev_end = ends[i]
            }
        }
    }
}' $ORIGINAL_GFF > $GFF_CLEAN

## call variants with minimap2 and iVar
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
        samtools view -@ 48 -bS ${OUTPUT_DIR}/${SAMPLE}.sam > ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
		samtools sort -@ 48 -o ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.bam
        samtools index ${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam  

        BAM="${OUTPUT_DIR}/${SAMPLE}.${POOLID}.sort.bam"

		## mpileup and ivar using gene annotations from the non-overlapping GFF
        samtools mpileup -A -d 0 -B -Q 0 -f $REF $BAM | \
		ivar variants -p ${OUTPUT_DIR}/${SAMPLE} -t 0 -m 1 -r $REF -g $GFF_CLEAN

		## keep only the non-empty tsv files
		## the tsv files contain headers so keep only the ones with more than one line
		if [ $(wc -l < "${OUTPUT_DIR}/${SAMPLE}.tsv") -le 1 ]; then
			rm -f "${OUTPUT_DIR}/${SAMPLE}.tsv"
		else
	    	echo "Variants found for ${SAMPLE}"
		fi
	done
fi

## keep the clean GFF and report the pool completed
echo "Clean GFF saved at: $GFF_CLEAN"
echo "Variant calling completed for pool ${POOLID}"