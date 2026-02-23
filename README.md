# mutatea
This is a framework comparing the mutational spectra of pathogen sequencing data between sources/cohorts.

`Inputs`:
1. metadata
2. either paired-end short reads or single reads (Illumina is default)

## crm: goal is to make this usable for reads of other sequencing types, remove limiting input (not just Illumina)
## crm: need to add an argument, minimap2 alignment depends on the type of sequencing read being run (can set default as Illumina seq since that's what we use)

3. reference fna(.gz) and gff(.gz)

`Outputs`:
## crm: output is unclear
1. tsv of mutations found in the data of each month 

2. optional: alignment files of the reads

3. optional: processed metadata files

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
```
git clone https://github.com/tiszalab/mutatea.git
```

2. Confirm installation of Python and required dependencies (pandas, biopython, openpyxl)

3. Activate conda environment with required packages

4. cd into mutatea folder and pip install `mutatea` package 
```
cd mutatea
```
```
pip install -e .
```

5. Confirm installation of function
```
mutatea -h
```

## crm: probably needs more steps

# Usage
```bash
mutatea -p <PATHOGEN_NAME> -m <WASTEWATER_METADATA_FILES> -pr <PAIRED_WASTEWATER_READS> -ref <REFERENCE_FILES>
```

# Required arguments
- `-p`, `--pathogen`: Pathogen to process, name should match the naming of the reads
- `-m`, `--wastewater_metadata`: Path to folder containing the wastewater metadata files (.xlsx)
- `-ref`, `--references`: Path to folder containing the reference files fna(.gz) and gff(.gz)
## crm: why do they have different file formats accepted for each read type? should be the same for both inputs
## crm: paired reads is not the accurate term
### Either
- `-pr`, `--paired_reads`: Path to folder containing paired wastewater reads (fastq)
- `-sr`, `--single_reads`: Path to folder containing single wastewater reads (fastq)
# crm: single reads can be fasta or fastq rn

# Optional arguments

## Data Configuration
- `-c`, `--clinical`: Path to folder containing the clinical metadata files (.xlsx) and fasta if parallel analysis is desired
- `-my`, `--monthly_only`: Only group wastewater samples by month, overrides default of being grouped by both month and month_region
- `-d`, `--dictionary`: Input custom mapping dictionary to map cities to any region (public health region, county, state, etc), overrides default of mapping cities to Texas public health regions
- `-g`, `--grouping`: Group samples by year, month, week, or day, overrides the default of samples grouped by month
- `-mw`, `--minimap_wastewater`: Overwrite minimap2 preset for the wastewater alignment (default is sr)
- `-mc`, `--minimap_clinical`: Overwrite minimap2 preset for the clinical alignment (default is asm10)

# crm: confirm that these are the default minimap2 presets I'm keeping (before publishing package)


## Output and Performance
- `-o`, `--output`: Path to desired output directory 
- `-f`, `--fast`: Override default for parallel workers, will run with all available cpus
- `-a`, `--all`: Keep all intermediate alignment files (otherwise deleted)
- `-l`, `--logger`: Export a detailed logger file
- `-s`, `--statistics`: See statistics on genome depth and coverage for each grouping

## Information
- `-tr`, `--timerange`: Print time range covered by the wastewater samples (clinical data should be time matched)
- `-v`, `--version`: Print the current version of mutatea



# Example
``` bash
mutatea -p H1N1 -m path/to/wastewater/metadata -pr /path/to/paired/wastewater/reads -ref path/to/ref/files -c path/to/clinical/files
```