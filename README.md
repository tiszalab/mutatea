# flu_cli
This is a framework comparing the mutational spectra of IAV (subtypes H1N1, H3N2, and H5N1) between sources/cohorts.

`Inputs`:
1. metadata
2. paired-end short reads
3. reference fna and gff

`Outputs`:
1. alignment files of the reads to the reference
2. tsv of mutations found in the reads

# Installation

1. Clone this github repository

2. Create the `conda` environment using the .yaml file, e.g.
`conda env create -f environment/flu_cli.yaml`

3. Activate the environment, i.e.
`conda activate flu_cli`

4. Use pip to install this command line tool, i.e.
`cd python_CLI`
`pip install .`

5. Usage
`flu_CLI -s <IAV_SUBTYPE> -m <WASTEWATER_METADATA_FILES> -r <WASTEWATER_READS> -ref <REFERENCE_FILES>'

# Required arguments
- '-s, --subtype': Influenza A subtype (H1N1, H3N2, H5N1)
- '-m, --wastewater_metadata': Path to folder containing the wastewater metadata files 
- '-r, --wastewater_reads': Path to folder containing the paired wastewater reads (fastq(.gz))
- '-ref, --reference_files': Path to folder containing the reference files fna(.gz) and gff(.gz)

optional arguments
- '-o, --output_dir': Path to output directory 
- '-c, --clinical_files': Path to folder containing the clinical metadata files and fasta
- `-my, --monthly_only': Only group wastewater samples by month (overrides default of being grouped by both month and month_region)
- '-t, --time_range': View time range of wastewater samples

# crm: file paths need to be cleaned here
# Example
`flu_CLI -s H1N1 -m /wastewater_metadata -r /path/to/wastewater/reads -ref /ref_files`