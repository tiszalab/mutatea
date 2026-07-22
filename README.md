<p align="center">
  <img src="mutatea_logo.png" width="300" alt="mutatea logo"/>
</p>

# mutatea
A framework for comparing the mutational spectra of pathogen sequencing data across sources and cohorts. Aligns  wastewater (and optionally clinical) reads to a reference genome, groups samples by time and/or region, and annotates all detected variants with coding effects using [varmint](https://github.com/tiszalab/varmint).

## Inputs
1. **Wastewater metadata** — one or more `.xlsx` files, each requires the following columns: `SampleID`, `City`, `Date` (YYYY-MM-DD, YYYY-MM, or YYYY)
2. **Wastewater reads** — accepted pre-aligned BAM files, paired reads, or single reads

# crm: clean up
paired-end or single-end reads (fastq/fasta, optionally gzipped) file names must contain the pathogen name


the pre-aligned BAM file does not require the pathogen name, but must have been aligned to the inputted reference genome

   - **Single read**: reads can have any of the following patterns — `<pathogen>.fasta`, `<pathogen>.fastq`, `<pathogen>.fastq.gz` 
     - Example: `sampleID.H1N1.fastq.gz`
   - **Paired reads**: R1/R2 files can be any of the following patterns  — `<pathogen>.R1.fastq.gz`, `<pathogen>.R1.fasta`, `<pathogen>_1.fastq`, `<pathogen>_1.fastq.gz`
     - Example: `sampleID.sars_cov2_R1.fastq`
   - **Pre-aligned BAM files**: The BAM file must have been aligned against the same reference genome, will only be returned if it contains reads found to have been aligned to contigs of the inputted reference genome
     - Example: `sampleID.sort.bam`
3. **Reference genome** — a folder containing one `.fna`/`.fna.gz` and one `.gff`/`.gff.gz` file
4. **Clinical sequences** *(optional)* — a folder containing `.fasta` files named by accession, plus a `.csv` metadata file with columns: `Accession`, `Collection_Date`

## Outputs
| Path | Description |
|------|-------------|
| `tsv_output/wastewater` | Per-group variant TSVs for wastewater sequences (time; time+region) |
| `tsv_output/clinical` | Per-group variant TSVs for clinical sequences |
| `alignment_files` | Merged BAMs per time group (and region) |
| `metadata_files` | Processed wastewater and clinical metadata CSVs |
| `<pathogen>_mutatea.log` | Detailed run log |
| `statistics` *(optional)* | samtools stats output per group |

### Output TSV columns
`contig`, `pos`, `var_type` (SNV/INS/DEL), `allele_type`, `ref_seq`, `alt_seq`, `depth`, `allele_count`, `AF`, `allele_avgq`, `allele_avgmq`, `ref_fwd`, `ref_rev`, `alt_fwd`, `alt_rev`, `strand_bias_p`, `ts_tv`, `VCF_PASS`, `is_coding`, `gene`, `transcript_id`, `strand`, `codon_ref`, `codon_alt`, `aa_ref`, `aa_alt`, `codon_index`, `codon_pos`, `effect`

# Installation
## Option A — conda + pip (recommended)
`mutatea.yaml` creates the conda environment (Python, minimap2, samtools, and Python dependencies); `pip install -e .` then installs the `mutatea` command via `pyproject.toml`.
```bash
git clone https://github.com/tiszalab/mutatea.git
cd mutatea
conda env create -f mutatea.yaml
conda activate mutatea_env
pip install -e .
```

## Option B — pip only
Installs the `mutatea` command and all Python dependencies via `pyproject.toml`.
Packages minimap2 and samtools must be available on your PATH separately (e.g. via `conda install -c bioconda minimap2 samtools`).
```bash
git clone https://github.com/tiszalab/mutatea.git
cd mutatea
pip install -e .
```

## Confirm installation
```bash
mutatea -h
```

# Usage
```bash
mutatea -p <PATHOGEN> -m <METADATA_DIR> -ref <REFERENCE_DIR> -d <DICTIONARY_JSON> -pr <PAIRED_READS_DIR>
```

# Required Arguments
- `-p`, `--pathogen`: Pathogen name — must match the naming convention used in the read files
- `-m`, `--wastewater_metadata`: Path to folder containing wastewater metadata files (`.xlsx`)
- `-ref`, `--references`: Path to folder containing reference `.fna`(`.gz`) and `.gff`(`.gz`) files
- `-d`, `--dictionary`: Path to a JSON file mapping cities to regions

One of the following read inputs is required:
- `-pr`, `--paired_reads`: Path to folder containing paired-end wastewater reads
- `-sr`, `--single_reads`: Path to folder containing single-end wastewater reads
- `-b`, `--bams`: Path to folder containing pre-aligned wastewater BAM files (must have been aligned to contigs of the reference genome you're running with)

# Optional Arguments

## Data Configuration
- `-c`, `--clinical`: Path to folder containing clinical fasta files and metadata CSV for parallel analysis
- `-ty`, `--time_only`: Group wastewater samples by time only, skipping time+region grouping
- `-g`, `--grouping`: Time grouping — `year`, `month`, `week`, or `day` (default: `month`)
- `-mw`, `--minimap_wastewater`: minimap2 preset for wastewater alignment (default: `sr`)
- `-mc`, `--minimap_clinical`: minimap2 preset for clinical alignment (default: `asm10`)
- `-q`, `--mapq`: Minimum mapping quality score for read filtering (default: `0`, no filtering)

## Output and Performance
- `-o`, `--output`: Path to output directory (default: current directory)
- `-f`, `--fast`: Use all available CPUs for parallel processing (recommended)
- `-a`, `--all`: Keep all intermediate alignment files (group-level BAMs are deleted by default after merging)
- `-s`, `--statistics`: Output per-group genome depth and coverage statistics

## Information
- `-tr`, `--timerange`: Print the date range covered by the wastewater samples
- `-v`, `--version`: Print the current version of mutatea

# Example
```bash
mutatea 
  -p H1N1 \
  -m path/to/wastewater/metadata \
  -ref path/to/reference/fasta/and/gff \
  -d path/to/json/dictionary \
  -pr path/to/paired/wastewater/reads \
  -c path/to/clinical/files \
  -q 20 \
  -f
```

# Dependencies
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://www.htslib.org/)
- [varmint](https://github.com/tiszalab/varmint) — variant calling and coding effect annotation
- pandas, biopython, pysam, openpyxl (installed automatically)

# License
MIT