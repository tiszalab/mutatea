# mutatea
A framework for comparing the mutational spectra of pathogen sequencing data across sources and cohorts (e.g. wastewater vs. clinical).

## Inputs
1. Wastewater metadata (.xlsx)
2. Paired-end or single-end wastewater reads (fasta/fastq)
3. Reference genome: fna(.gz) and gff(.gz)
4. Optional: clinical fasta files for parallel analysis

## Outputs
1. TSV files describing mutations found in each grouped sample set
2. Optional: merged alignment files (BAM) per time group and region
3. Optional: genome coverage and depth statistics
4. Optional: processed metadata files

# Installation

## Prerequisites
Install the following tools via conda (recommended):
```bash
conda install -c bioconda minimap2 samtools
```

## Steps
1. Clone this repository:
```bash
git clone https://github.com/tiszalab/mutatea.git
```

2. Install the Python package and its dependencies (pandas, biopython, pysam, openpyxl):
```bash
cd mutatea
pip install -e .
```

3. Confirm installation:
```bash
mutatea -h
```

# Usage
```bash
mutatea -p <PATHOGEN> -m <METADATA_DIR> -pr <PAIRED_READS_DIR> -ref <REFERENCE_DIR>
```

# Required Arguments
- `-p`, `--pathogen`: Pathogen name — must match the naming convention used in the read files
- `-m`, `--wastewater_metadata`: Path to folder containing wastewater metadata files (.xlsx)
- `-ref`, `--references`: Path to folder containing reference fna(.gz) and gff(.gz) files

One of the following read inputs is required:
- `-pr`, `--paired_reads`: Path to folder containing paired-end wastewater reads
- `-sr`, `--single_reads`: Path to folder containing single-end wastewater reads

# Optional Arguments

## Data Configuration
- `-c`, `--clinical`: Path to folder containing clinical fasta files for parallel analysis
- `-ty`, `--time_only`: Group wastewater samples by time only, skipping time+region grouping
- `-d`, `--dictionary`: Custom mapping dictionary to assign cities to regions (default: Texas public health regions)
- `-g`, `--grouping`: Time grouping resolution — `year`, `month`, `week`, or `day` (default: `month`)
- `-mw`, `--minimap_wastewater`: minimap2 preset for wastewater alignment (default: `sr`)
- `-mc`, `--minimap_clinical`: minimap2 preset for clinical alignment (default: `asm10`)
- `-q`, `--mapq`: Minimum mapping quality score for read filtering (default: `0`, no filtering)

## Output and Performance
- `-o`, `--output`: Path to output directory
- `-f`, `--fast`: Use all available CPUs for parallel processing
- `-a`, `--all`: Keep all intermediate alignment files (pool-level BAMs are deleted by default after merging)
- `-l`, `--logger`: Write a detailed log file
- `-s`, `--statistics`: Output per-group genome depth and coverage statistics

## Information
- `-tr`, `--timerange`: Print the time range covered by the wastewater samples
- `-v`, `--version`: Print the current version of mutatea

# Example
```bash
mutatea -p H1N1 -m path/to/wastewater/metadata -pr /path/to/paired/wastewater/reads -ref path/to/ref/files -c path/to/clinical/files
```