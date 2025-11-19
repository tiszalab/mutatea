#!/bin/bash

# Script to download Influenza A H1N1 complete genomes from humans in Texas
# collected between 2022-05-04 and 2025-11-05 using Entrez Direct.

# require flu subtype
SUBTYPE=$1   

# User email for NCBI usage policy
EMAIL="u255582@bcm.edu"
export NCBI_EMAIL="$EMAIL"

# Output file
OUTFILE="${SUBTYPE}_clinical.fasta"

# EDirect search and fetch
esearch -db nuccore \
  -query "Influenza A virus[Organism] AND txid11320[Organism] \
  AND complete genome[Title] \
  AND Homo sapiens[Host] AND txid9606[Host] \
  AND 2022/05/04:2025/11/05[CollectionDate] \
  AND ${SUBTYPE}[Genotype] \
  AND USA: TX[Country]" \
  | efetch -format fasta > "${OUTFILE}"

echo "Download complete! Sequences saved to ${OUTFILE}"
