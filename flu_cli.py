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
import tempfile

# load in functions from metadata_funcs
# crm: make sure to update names of functions being imported as they change in metadata_funcs.py
try:
    from .metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_merge_key_lists
except:
    from metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_merge_key_lists

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
    # set default file path for output to the current working directory of the user
    output_path_default = os.getcwd()

    # create parser and describe function of script
    parser = argparse.ArgumentParser(description="Process and align wastewater and clinical influenza A reads for mutation analysis.")

    ## required arguments
    # argument for IAV subtype
    parser.add_argument("-s", "--subtype", type=str, required=True, choices=["H1N1", "H3N2", "H5N1"], help="Influenza subtype to process, options are H1N1, H3N2, and H5N1")

    # argument for file path to folder containing wastewater metadata files
    parser.add_argument("-m", "--wastewater_metadata", type=str, required=True, help="Path to folder containing wastewater metadata files")

    # crm: currently set reads to DNR, will update after I can start accessing the reads
    # argument for file path to folders containing wastewater reads
    parser.add_argument("-r", "--wastewater_reads", type=str, required=True, help="Path to the folders containing the wastewater reads")

    # argument for file path to folder containing reference files
    parser.add_argument("-ref", "--reference_files", type=str, required=True, help="Path to folder containing the reference fasta(.gz) and gff(.gz) files")

    ## optional arguments
    # argument for file path to folder containing clinical files
    parser.add_argument("-c", "--clinical_files", type=str, help="Path to folder containing the clinical fasta and csv files")

    # argument for nondefault output directory
    parser.add_argument("-o", "--output_dir", type=str, default=output_path_default, help="Path to chosen output directory")

    # argument to only split wastewater metadata by month
    parser.add_argument("-my", "--month_only", action='store_true', help="Override default split of month and region, will only split the wastewater data by month")

    # argument to view time range covered by wastewater metadata
    parser.add_argument("-t", "--time_range", type=str, help="View time range covered by wastewater sample collection")

    # parse arguments
    args = parser.parse_args()

    ## boolean checks
    # check if clinical files are included
    include_clinical = bool(args.clinical_files)

    # check if region is included
    include_region = not args.month_only

    # initialize directories dictionary
    dirs = {}
    
    # create main output directory with subtype-specific subfolder
    dirs["output"] = os.path.join(args.output_dir, f"{args.subtype}_align")
    os.makedirs(dirs["output"], exist_ok=True)

    # define logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)

    # crm: not sure if I also need filehandler, but it would be put in this section fyi

    logger.addHandler(stream_handler)
    
    ############################## process metadata files ##############################
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
    
    # crm: messy, clean this up
    # optionally give the time range of the wastewater samples
    if args.time_range:
        logger.info(f"\nViewing the time range of the wastewater samples")
        get_date_range(metadata)

    # create directory for cleaned and merged metadata
    dirs["metadata_dir"] = os.path.join(dirs["output"], "metadata_files")
    os.makedirs(dirs["metadata_dir"], exist_ok=True)

    # export cleaned wastewater metadata
    logger.info(f"\nExporting the cleaned wastewater metadata to {dirs['metadata_dir']}")
    export_metadata(metadata, dirs["metadata_dir"])

    # load in clinical metadata
    if include_clinical:
        logger.info("\nLoading in clinical metadata")
        clinical_metadata = load_clinical_metadata(args.clinical_files)

        # export cleaned clinical metadata
        logger.info(f"\nExporting the cleaned clinical metadata to {dirs['metadata_dir']}")
        clinical_metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_clinical_{args.subtype}.csv"), index=False)
    
   
    # create directory for unzipped reference files
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    ## process reference files
    logger.info("\nProcessing reference files")
    reference_files = process_reference_file(args.reference_files, dirs["reference_dir"])

    ############################## process reads ##############################
    # find wastewater reads from pools
    logger.info("\nFinding wastewater reads from pools")
    wastewater_reads = find_wastewater_reads(args.wastewater_reads, args.subtype)
    
    # create file directory for alignment files
    dirs["alignment_dir"] = os.path.join(dirs["output"], "alignment_files")
    os.makedirs(dirs["alignment_dir"], exist_ok=True)

    # create subfolders in the alignment directory
    dirs["wastewater_dir"] = os.path.join(dirs["alignment_dir"], "wastewater")
    os.makedirs(dirs["wastewater_dir"], exist_ok=True)

    # create directory for pools
    dirs["pools"] = os.path.join(dirs["wastewater_dir"], "pools")
    os.makedirs(dirs["pools"], exist_ok=True)
    
    # align wastewater reads to reference genome
    logger.info("\nAligning wastewater reads to reference genome")
    align_wastewater_reads(wastewater_reads, dirs["reference_dir"], dirs["pools"])

    # create directory for wastewater lists
    dirs["wastewater_lists_dir"] = os.path.join(dirs["wastewater_dir"], "lists")
    os.makedirs(dirs["wastewater_lists_dir"], exist_ok=True)

    # create subfolders for wastewater lists
    if include_region == True:
        dirs["wastewater_list_month"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month")
        os.makedirs(dirs["wastewater_list_month"], exist_ok=True)
        dirs["wastewater_list_region"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month_region")
        os.makedirs(dirs["wastewater_list_region"], exist_ok=True)
    else:
        dirs["wastewater_list_month"] = dirs["wastewater_lists_dir"]
        os.makedirs(dirs["wastewater_list_month"], exist_ok=True)

    # crm: set to DNR for now
    # wastewater merged bams
    #dirs["merged_bams"] = os.path.join(dirs["wastewater_dir"], "merged_bams")
    #os.makedirs(dirs["merged_bams"], exist_ok=True)

    # create subfolder: wastewater bams merged by month
    dirs["merged_bams_month"] = os.path.join(dirs["merged_bams"], "merged_bams_month")        
    os.makedirs(dirs["merged_bams_month"], exist_ok=True)

    # crm: need to add function for merging bams by month

    # create subfolder: wastewater bams merged by month and region
    if include_region:
        dirs["merged_bams_month_region"] = os.path.join(dirs["merged_bams"], "merged_bams_month_region")
        os.makedirs(dirs["merged_bams_month_region"], exist_ok=True)

    # crm: need to add function for merging bams by month and region (if region included)
    
    # crm: set to DNR for now
    # create tsv_output folder to later catch tsv files
    #dirs["tsv_output"] = os.path.join(dirs["output"], "tsv_output")
    #os.makedirs(dirs["tsv_output"], exist_ok=True)

    # crm: need to add in varmint for wastewater

    

    # process clinical files
    if include_clinical:
        # create folder for clinical output
        dirs["clinical"] = os.path.join(dirs["alignment_dir"], "clinical")
        os.makedirs(dirs["clinical"], exist_ok=True)

        # load in clinical fasta
        logger.info("\nLoading in clinical FASTA")
        clinical_fasta = load_clinical_fasta(args.clinical_files)

        # create folder for the lists of accessions by month
        dirs["clinical_lists_month"] = os.path.join(dirs["clinical"], "clinical_lists_month")
        os.makedirs(dirs["clinical_lists_month"], exist_ok=True)

        # create monthly lists of accessions
        logger.info("\nCreating monthly lists of accessions")
        create_monthly_accession_lists(clinical_metadata, dirs["clinical_lists_month"])

        # crm: want to later remove these clinical fasta files, want to use tempfile to instead save them to a temp dir
        # create folder for the clinical fastas split by month
        dirs["clinical_fasta_month"] = os.path.join(dirs["clinical"], "clinical_fasta_month")
        os.makedirs(dirs["clinical_fasta_month"], exist_ok=True)

        # split clinical fasta by monthly lists
        logger.info("\nSplitting clinical FASTA by monthly lists")
        split_clinical_fasta_by_month(clinical_fasta, dirs["clinical_lists_month"], dirs["clinical_fasta_month"])

        # create folder for the clinical bam files that were merged by month
        # dirs["clinical_bam_month"] = os.path.join(dirs["clinical"], "clinical_bam_month")
        # os.makedirs(dirs["clinical_bam_month"], exist_ok=True)

        # crm: add function to align clinical reads to reference
        # crm: output should we sorted bam files

        # crm: varmint for clinical bam files
        # crm: output should be tsv files saved to tsv_output


    # CRM: would like to add in a "find existing reference files" function to metadata_funcs.py
        