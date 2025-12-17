#!/bin/bash
set -euo pipefail

# Directories
BAM_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/p1700"     # directory with all your sorted BAMs
REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference.fasta"
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/GCA_039113225.1_ASM3911322v1_genomic.gff"
OUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/ivar_variants"

mkdir -p "$OUT_DIR"

# Loop over BAMs
for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "Calling variants for sample $SAMPLE..."

    samtools mpileup -aa -A -d 0 -B -Q 0 -f "$REF" "$BAM" | \
        ivar variants \
            -p "$OUT_DIR/${SAMPLE}" \
            -t 0.03 \
            -m 10 \
            -q 20 \
            -r "$REF" \
            -g "$GFF"
done

echo "Variant calling complete. TSVs are in $OUT_DIR"
