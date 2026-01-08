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
from datetime import timedelta

# load in functions from metadata_funcs
# crm: make sure to update names of functions being imported as they change in metadata_funcs.py
try:
    from .metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_lists, merge_wastewater_bams, align_clinical_reads, varmint
except:
    from metadata_funcs import load_merge_metadata, add_month_year, add_region, ensure_sitecode_column, reorganize_metadata_columns, export_metadata, get_date_range, load_clinical_metadata, load_clinical_fasta, process_reference_file, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_lists, merge_wastewater_bams, align_clinical_reads, varmint

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
    parser.add_argument("-t", "--time_range", action='store_true', help="View time range covered by wastewater sample collection")

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
    metadata = load_merge_metadata(args.wastewater_metadata)
    if metadata.empty:
        logger.error("No metadata files found in the specified directory.")
        sys.exit(1)

    # add month_year column to metadata
    metadata = add_month_year(metadata)

    # add region column to the metadata if user didn't specify month only
    if not args.month_only:
        # add region column to metadata
        metadata = add_region(metadata)
        no_region = False
    else:
        no_region = True

    # confirm SiteCode column is in metadata
    metadata = ensure_sitecode_column(metadata)
    
    # reorganize metadata columns with consideration of region potentially not being included
    metadata = reorganize_metadata_columns(metadata, no_region=no_region)
    
    # optionally give time range of the wastewater samples
    if args.time_range:
        earliest, latest = get_date_range(metadata)
        logger.info(f"Date range covered by wastewater data: {earliest.strftime('%Y-%m-%d')} to {latest.strftime('%Y-%m-%d')}")

    # create directory for processed and merged metadata
    dirs["metadata_dir"] = os.path.join(dirs["output"], "metadata_files")
    os.makedirs(dirs["metadata_dir"], exist_ok=True)

    # export processed wastewater metadata
    logger.info(f"Exporting the processed metadata to {dirs['metadata_dir']}\n")
    export_metadata(metadata, dirs["metadata_dir"])

    # load in clinical metadata
    if include_clinical:
        clinical_metadata = load_clinical_metadata(args.clinical_files)

        # export processed clinical metadata
        clinical_metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_clinical_{args.subtype}.csv"), index=False)
    
   
    # create directory for unzipped reference files
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    ## process reference files
    reference_files = process_reference_file(args.reference_files, dirs["reference_dir"])

    ############################## process reads ##############################
    # find wastewater reads from pools
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
    logger.info("\nAligning wastewater reads to given reference genome\n")
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

    # create list of monthly accessions for downstream merge key creation
    if include_region == True:
        logger.info("\nMerging wastewater alignment files by month and public health region")
    else:
        logger.info("\nMerging wastewater alignment files by month")

    create_wastewater_bam_lists(metadata, dirs["pools"], dirs["wastewater_list_month"], dirs.get("wastewater_list_region"), include_region)

    # wastewater merged bams
    dirs["merged_bams"] = os.path.join(dirs["wastewater_dir"], "merged_bams")
    os.makedirs(dirs["merged_bams"], exist_ok=True)

    # create subfolder: wastewater bams merged by month
    dirs["merged_bams_month"] = os.path.join(dirs["merged_bams"], "merged_bams_month")        
    os.makedirs(dirs["merged_bams_month"], exist_ok=True)

    # merge wastewater bams by month
    merge_wastewater_bams(dirs["wastewater_list_month"], dirs["merged_bams_month"])

    # create subfolder: wastewater bams merged by month and region
    if include_region:
        dirs["merged_bams_month_region"] = os.path.join(dirs["merged_bams"], "merged_bams_month_region")
        os.makedirs(dirs["merged_bams_month_region"], exist_ok=True)
        
        # merge wastewater bams by month and region
        merge_wastewater_bams(dirs["wastewater_list_region"], dirs["merged_bams_month_region"])
    
    # create tsv_output folder to later catch tsv files
    dirs["tsv_output"] = os.path.join(dirs["output"], "tsv_output")
    os.makedirs(dirs["tsv_output"], exist_ok=True)

    # create subfolders if clinical included
    if include_clinical:
        # split output tsv by source
        dirs["tsv_wastewater"] = os.path.join(dirs["tsv_output"], "wastewater")
        os.makedirs(dirs["tsv_wastewater"], exist_ok=True)

        dirs["tsv_clinical"] = os.path.join(dirs["tsv_output"], "clinical")
        os.makedirs(dirs["tsv_clinical"], exist_ok=True)
        
        # split clinical output by grouping method
        if include_region:
            dirs["tsv_month"] = os.path.join(dirs["tsv_wastewater"], "month")
            os.makedirs(dirs["tsv_month"], exist_ok=True)
            dirs["tsv_month_region"] = os.path.join(dirs["tsv_wastewater"], "month_region")
            os.makedirs(dirs["tsv_month_region"], exist_ok=True)

    # or just split clinical output by grouping method
    else:
        dirs["tsv_month"] = os.path.join(dirs["tsv_output"], "month")
        os.makedirs(dirs["tsv_month"], exist_ok=True)

        dirs["tsv_month_region"] = os.path.join(dirs["tsv_output"], "month_region")
        os.makedirs(dirs["tsv_month_region"], exist_ok=True)
            
    # crm: still in test
    # varmint for wastewater
    logger.info("\nAnnotating coding effects of mutations with varmint")

    # crm: this is too messy, please find a better way to do this
    if include_region:
        varmint(dirs["merged_bams_month"], dirs["reference_dir"], dirs["tsv_month"])
        varmint(dirs["merged_bams_month_region"], dirs["reference_dir"], dirs["tsv_month_region"])
    else:
        varmint(dirs["merged_bams_month"], dirs["reference_dir"], dirs["tsv_output"])


    














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
        # crm: I want to delete this logger line, which I could go forward with if I am saving the monthly lists to tempdir
        #logger.info("\nCreating monthly lists of accessions")
        create_monthly_accession_lists(clinical_metadata, dirs["clinical_lists_month"])

        # crm: want to later remove these clinical fasta files, want to use tempfile to instead save them to a temp dir
        # create folder for the clinical fastas split by month
        dirs["clinical_fasta_month"] = os.path.join(dirs["clinical"], "clinical_fasta_month")
        os.makedirs(dirs["clinical_fasta_month"], exist_ok=True)

        # split clinical fasta by monthly lists
        logger.info("\nSplitting clinical FASTA by month")
        split_clinical_fasta_by_month(clinical_fasta, dirs["clinical_lists_month"], dirs["clinical_fasta_month"])

        # create folder for the clinical bam files that were merged by month
        dirs["clinical_bam_month"] = os.path.join(dirs["clinical"], "clinical_bam_month")
        os.makedirs(dirs["clinical_bam_month"], exist_ok=True)

        # align clinical reads to reference
        logger.info("\nAligning clinical reads to the reference genome")
        align_clinical_reads(dirs["clinical_fasta_month"], dirs["clinical_bam_month"], dirs["reference_dir"])

        # varmint for clinical
        logger.info("\nAnnotating coding effects of mutations with varmint")
        varmint(dirs["clinical_bam_month"], dirs["reference_dir"], dirs["tsv_clinical"])


    # crm: pretty sure I needed more things at the end than this
    # print run time
    cli_end_time = time.perf_counter()
    time_taken = round((cli_end_time - cli_start_time), 2) 
    # crm: can make this prettier later
    logger.info(f"Run time: {timedelta(seconds=time_taken)}")