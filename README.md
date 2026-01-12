# crm: virmuth
This is a framework comparing the mutational spectra of virome sequencing data between sources/cohorts.

`Inputs`:
1. metadata
2. paired-end short reads
    
    crm: could be one single read too, e.g. COVID
3. reference fna and gff

`Outputs`:
1. alignment files of the reads to the reference
2. tsv of mutations found in the reads

# Installation
## Prerequisites
Must install the following tools:
- **minimap2** (for read alignment)
- **samtools** (for BAM file processing)
- **varmint** (for annotating coding effects using GFF CDS features)

Install via conda (recommended):
```bash
conda install -c bioconda minimap2 samtools
```

Quick install varmint with `pip`

*This will not install dependencies, it will be installed into whichever environment you are in*
```bash
conda install pip
cd /path/to/varmint
pip install .
```

## Steps
1. Clone this github repository

2. Install the Python package and dependencies (pandas, biopython, openpyxl)

3. Install varmint and dependencies

4. Activate conda environment with required packages

# crm: probably needs more steps

# Usage
```bash
flu_CLI -s <IAV_SUBTYPE> -m <WASTEWATER_METADATA_FILES> -r <WASTEWATER_READS> -ref <REFERENCE_FILES>
```

# Required arguments
- `-s`, `--subtype`: Influenza A subtype (H1N1, H3N2, H5N1) #crm: need to adjust
- `-m`, `--wastewater_metadata`: Path to folder containing the wastewater metadata files 
- `-r`, `--wastewater_reads`: Path to folder containing the paired wastewater reads (fastq(.gz))
- `-ref`, `--reference_files`: Path to folder containing the reference files fna(.gz) and gff(.gz)

# Optional arguments
- `-o`, `--output_dir`: Path to output directory 
- `-c`, `--clinical_files`: Path to folder containing the clinical metadata files and fasta
- `-my`, `--monthly_only`: Only group wastewater samples by month (overrides default of being grouped by both month and month_region)
- `-t`, `--time_range`: View time range of wastewater samples
- `-v`, `--version`: View the current version of crm: flu cli

# Example
`flu_CLI -s H1N1 -m path/to/wastewater/metadata -r /path/to/wastewater/reads -ref path/to/ref/files`