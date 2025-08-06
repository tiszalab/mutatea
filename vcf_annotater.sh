#!/bin/bash

## a script for adding gene annotations to VCF files using bcftools
## requires a VCF file and a GFF3 file

VCF=$1

REF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/H1N1_reference.fasta"
GFF="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/GCA_039113225.1_ASM3911322v1_genomic.gff"


bcftools csq -f $REF -g $GFF $VCF -Ob -o annotated.bcf
bcftools view annotated.bcf -Oz -o annotated.vcf.gz
bcftools index annotated.vcf.gz