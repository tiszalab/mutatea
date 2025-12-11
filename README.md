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
`cd flu_cli`
`pip install .`
