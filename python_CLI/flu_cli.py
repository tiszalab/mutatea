#!/usr/bin/env python

###################### SETUP ######################
# load modules
import argparse
import sys, os
import logging
import gzip
import shutil
from Bio import SeqIO
from pathlib import Path
import time

# load in functions from metadata_funcs
# crm: make sure to update names of functions being imported as they change in metadata_funcs.py
try:
    from .metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, create_monthly_accession_lists, split_clinical_fasta_by_month, create_output_directories
except:
    from metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, create_monthly_accession_lists, split_clinical_fasta_by_month, create_output_directories

# convert string to boolean for argparse
def str2bool(x):
    if isinstance(x, bool):
       return x
    if x.lower() in ("yes", "y"):
        return True
    elif x.lower() in ('no', 'n'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# entry point function for the CLI
def flu_cli():
    # start timer
    cli_start_time = time.perf_counter()
    # set default file path for output to be the folder where this python script is kept
    output_path_default = os.path.dirname(os.path.abspath(__file__))

    # create parser and describe function of script
    parser = argparse.ArgumentParser(description="Process and align wastewater and clinical influenza A reads for mutation analysis.")

    ## required arguments
    # argument for IAV subtype
    parser.add_argument("-s", "--subtype", type=str, required=True, choices=["H1N1", "H3N2", "H5N1"], help="Influenza subtype to process, options are H1N1, H3N2, and H5N1")

    # argument for file path to folder containing wastewater metadata files
    parser.add_argument("-m", "--wastewater_metadata", type=str, required=True, help="Path to folder containing wastewater metadata files")

    # argument for file path to folders containing wastewater reads
    # crm: add filter in help for formatting (e.g. subtype.R1.fastq)
    # parser.add_argument("-r", "--wastewater_reads", type=str, required=True, help="Path to the folders containing the wastewater reads")

    # argument for file path to folder containing reference files
    parser.add_argument("-ref", "--reference_files", type=str, required=True, help="Path to folder containing the reference fasta(.gz) and gff(.gz) files")

    ## optional arguments
    # argument for file path to folder containing clinical files
    parser.add_argument("-c", "--clinical_files", type=str, help="Path to folder containing the clinical fasta and csv files")

    # argument for nondefault output directory
    parser.add_argument("-o", "--output_dir", type=str, default=output_path_default, help="Path to chosen output directory")

    # argument to only split wastewater metadata by month
    parser.add_argument("-my", "--month_only", type=str, help="Override default split of month and region, only splits wastewater metadata by month")

    # argument to view time range covered by wastewater metadata
    parser.add_argument("-t", "--time_range", type=str, help="View time range covered by wastewater sample collection")

    # parse arguments
    args = parser.parse_args()

    # check if clinical files are included
    include_clinical = bool(args.clinical_files)

    ## make directories
    dirs = create_output_directories(args.output_dir, include_clinical)
    