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
    from .metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, create_output_directories
except:
    from metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, create_output_directories

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
    parser.add_argument("-my", "--month_only", type=str, help="Override default split of month and region, will only split the wastewater data by month")

    # argument to view time range covered by wastewater metadata
    parser.add_argument("-t", "--time_range", type=str, help="View time range covered by wastewater sample collection")

    # parse arguments
    args = parser.parse_args()

    # check if clinical files are included
    include_clinical = bool(args.clinical_files)

    # check if region is included
    include_region = not args.month_only

    # make directories
    dirs = create_output_directories(args.output_dir, include_region, include_clinical)

    # define logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)

    # crm: not sure if I need filehandler, but it would be put in this section fyi

    logger.addHandler(stream_handler)
    
    # process wastewater metadata
    logger.info(f"\nProcessing wastewater metadata from: {args.wastewater_metadata}")
    metadata = load_merge_metadata(args.wastewater_metadata)
    if metadata.empty:
        logger.error("No metadata files found in the specified directory.")
        sys.exit(1)
    
    # add month_year column to metadata
    logger.info("\nAdding month_year column to metadata")
    metadata = add_month_year(metadata)

    # add region column to the metadata if user didn't specify month only
    if not args.month_only:
        # add region column to metadata
        logger.info("\nAdding region column to metadata")
        metadata = add_region(metadata)
        no_region = False
    else:
        no_region = True

    # confirm SiteCode column is in metadata
    # crm: do I really need to say this to the user?
    logger.info("\nAdding SiteCode column to metadata")
    metadata = ensure_sitecode_column(metadata)
    
    # reorganize metadata columns with consideration of region potentially not being included
    logger.info("\nReorganizing metadata columns")
    metadata = reorganize_metadata_columns(metadata, no_region=no_region)
    
    # optionally give the time range of the wastewater samples
    if args.time_range:
        logger.info(f"\nViewing the time range of the wastewater samples")
        get_date_range(metadata)

    # export cleaned wastewater metadata
    logger.info(f"\nExporting the cleaned wastewater metadata to {dirs['metadata_dir']}")
    export_metadata(metadata, dirs["metadata_dir"])

    # load in clinical metadata
    if include_clinical:
        logger.info("\nLoading in clinical metadata")
        clinical_metadata = load_clinical_metadata(args.clinical_files)

    # crm: maybe move lower
    # load in clinical fasta
    if include_clinical:
        logger.info("\nLoading in clinical FASTA")
        clinical_fasta = load_clinical_fasta(args.clinical_files)

    # create monthly lists of accessions
    if include_clinical:
        logger.info("\nCreating monthly lists of accessions")
        create_monthly_accession_lists(clinical_metadata, dirs["clinical_lists_month"])
       
    # split clinical fasta by monthly lists
    if include_clinical:
        logger.info("\nSplitting clinical FASTA by monthly lists")
        split_clinical_fasta_by_month(args.clinical_files, dirs["clinical_lists_month"], dirs["clinical_fasta_month"])

    # CRM: find existing reference files

    # process reference files
    logger.info("\nProcessing reference files")
    reference_files = process_reference_file(args.reference_files, dirs["reference_dir"])