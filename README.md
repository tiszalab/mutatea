# mutatea
A framework for comparing the mutational spectra of pathogen sequencing data across sources and cohorts. mutatea aligns wastewater (and optionally clinical) reads to a reference genome, groups samples by time and/or region, and annotates all detected variants with coding effects using [varmint](https://github.com/tiszalab/varmint).

## Inputs
1. **Wastewater metadata** — one or more `.xlsx` files, each requiring columns: `SampleID`, `Date` (YYYY-MM-DD, YYYY-MM, or YYYY), `City`
2. **Wastewater reads** — paired-end or single-end reads (fastq/fasta, optionally gzipped); file names must contain the pathogen name (e.g. `sample.H1N1.R1.fastq`)
3. **Reference genome** — a folder containing one `.fna`/`.fna.gz` and one `.gff`/`.gff.gz` file
4. **Clinical sequences** *(optional)* — a folder containing `.fasta` files named by accession, plus a `.csv` metadata file with columns: `Accession`, `Collection_Date`

## Outputs
| Path | Description |
|------|-------------|
| `tsv_output/wastewater/` | Per-group variant TSVs (time; time+region) |
| `tsv_output/clinical/` | Per-group variant TSVs for clinical sequences |
| `alignment_files/` | Merged BAMs per time group (and region) |
| `metadata_files/` | Processed wastewater and clinical metadata CSVs |
| `statistics/` *(optional)* | samtools stats output per group |
| `*_mutatea.log` *(optional)* | Detailed run log |

### Output TSV columns
`contig`, `pos`, `var_type` (SNV/INS/DEL), `allele_type`, `ref_seq`, `alt_seq`, `depth`, `allele_count`, `allele_avgq`, `allele_avgmq`, `strand_bias_p`, `VCF_PASS`, `is_coding`, `gene`, `transcript_id`, `strand`, `codon_ref`, `codon_alt`, `aa_ref`, `aa_alt`, `codon_index`, `codon_pos`, `effect`

# Installation

## Option A — conda environment (recommended)
```bash
git clone https://github.com/tiszalab/mutatea.git
cd mutatea
conda env create -f mutatea.yaml
conda activate mutatea
pip install -e .
```

## Option B — pip only
```bash
git clone https://github.com/tiszalab/mutatea.git
cd mutatea
pip install -e .
```
> minimap2 and samtools must be available on your PATH (e.g. via `conda install -c bioconda minimap2 samtools`).

## Confirm installation
```bash
mutatea -h
```

# Usage
```bash
mutatea -p <PATHOGEN> -m <METADATA_DIR> -pr <PAIRED_READS_DIR> -ref <REFERENCE_DIR>
```

# Required Arguments
- `-p`, `--pathogen`: Pathogen name — must match the naming convention used in the read files
- `-m`, `--wastewater_metadata`: Path to folder containing wastewater metadata files (`.xlsx`)
- `-ref`, `--references`: Path to folder containing reference `.fna`(`.gz`) and `.gff`(`.gz`) files

One of the following read inputs is required:
- `-pr`, `--paired_reads`: Path to folder containing paired-end wastewater reads
- `-sr`, `--single_reads`: Path to folder containing single-end wastewater reads

# Optional Arguments

## Data Configuration
- `-c`, `--clinical`: Path to folder containing clinical fasta files and metadata CSV for parallel analysis
- `-ty`, `--time_only`: Group wastewater samples by time only, skipping time+region grouping
- `-d`, `--dictionary`: Path to a JSON file mapping city names to regions (default: Texas public health regions)
- `-g`, `--grouping`: Time grouping resolution — `year`, `month`, `week`, or `day` (default: `month`)
- `-mw`, `--minimap_wastewater`: minimap2 preset for wastewater alignment (default: `sr`)
- `-mc`, `--minimap_clinical`: minimap2 preset for clinical alignment (default: `asm10`)
- `-q`, `--mapq`: Minimum mapping quality score for read filtering (default: `0`, no filtering)

## Output and Performance
- `-o`, `--output`: Path to output directory (default: current directory)
- `-f`, `--fast`: Use all available CPUs for parallel processing
- `-a`, `--all`: Keep all intermediate alignment files (pool-level BAMs are deleted by default after merging)
- `-l`, `--logger`: Write a detailed log file to the output directory
- `-s`, `--statistics`: Output per-group genome depth and coverage statistics

## Information
- `-tr`, `--timerange`: Print the date range covered by the wastewater samples
- `-v`, `--version`: Print the current version of mutatea

# Example
```bash
mutatea -p H1N1 \
  -m path/to/wastewater/metadata \
  -pr path/to/paired/wastewater/reads \
  -ref path/to/reference/files \
  -c path/to/clinical/files \
  -q 20 -f -l
```

# Dependencies
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://www.htslib.org/)
- [varmint](https://github.com/tiszalab/varmint) — variant calling and coding effect annotation
- pandas, biopython, pysam, openpyxl (installed automatically)

# License
MIT
