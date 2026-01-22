# crm: virmuth
This is a framework comparing the mutational spectra of virome sequencing data between sources/cohorts.

`Inputs`:
1. metadata
2. either paired-end short reads or single reads
## crm: we don't use paired-end short reads, ONT doesn't do that
## crm: need to add an argument, minimap2 alignment depends on the type of sequencing read being run (can set default as ONT since that's what we use)

3. reference fna and gff

`Outputs`:
## crm: output is unclear
1. tsv of mutations found in the data of each month 
2. optional: alignment files of the reads

# Installation
## Prerequisites
Must install the following tools:
- **minimap2** (for read alignment)
- **samtools** (for BAM file processing)

Install via conda (recommended):
```bash
conda install -c bioconda minimap2 samtools
```

## Steps
1. Clone this github repository

2. Install the Python package and dependencies (pandas, biopython, openpyxl)

3. Activate conda environment with required packages

## 4. crm: pip install package and run as function

## crm: probably needs more steps

# Usage
```bash
flu_CLI -s <VIRUS_NAME> -m <WASTEWATER_METADATA_FILES> -pr <PAIRED_WASTEWATER_READS> -ref <REFERENCE_FILES>
```

# Required arguments
## crm: need to adjust, it's not just IAV subtype

- `-s`, `--subtype`: Influenza A subtype (H1N1, H3N2, H5N1) ????
- `-m`, `--wastewater_metadata`: Path to folder containing the wastewater metadata files (.xlsx)
- `-ref`, `--reference_files`: Path to folder containing the reference files fna(.gz) and gff(.gz)
## crm: why do they have different file formats accepted for each read type? should be the same for both inputs
## crm: paired reads is not the accurate term
### Either
- `-pr`, `--paired_reads`: Path to folder containing paired wastewater reads (fastq)
- `-sr`, `--single_reads`: Path to folder containing single wastewater reads (fasta)

# Optional arguments
- `-o`, `--output_dir`: Path to desired output directory 
- `-c`, `--clinical_files`: Path to folder containing the clinical metadata files (.xlsx) and fasta if parallel analysis is desired
- `a`, `--all`: Keep all intermediate alignment files (otherwise deleted)
- `-my`, `--monthly_only`: Only group wastewater samples by month (overrides default of being grouped by both month and month_region)
- `-t`, `--time_range`: Print time range of wastewater samples
- `-v`, `--version`: Print the current version of ????
## crm: clean up description of --version above

# Example
``` bash
flu_CLI -s H1N1 -m path/to/wastewater/metadata -pr /path/to/paired/wastewater/reads -ref path/to/ref/files -c path/to/clinical/files
```