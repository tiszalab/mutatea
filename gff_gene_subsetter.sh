#!/bin/bash

## a script for breaking down a GFF file to only include one gene per file
## this is needed for using ivar to get tsv files for each gene

GFF_ALL="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/GCA_039113225.1_ASM3911322v1_genomic.gff"
OUTPUT_DIR="/gpfs1/projects/Tisza_Lab/crm_flu_mutatome/H1N1_align/reference/gff/gene_gff"
