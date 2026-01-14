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
from datetime import timedelta

# load in functions from cli funcs
# crm: make sure to update names of functions being imported as they change in cli_funcs.py
try:
    from .cli_funcs import process_reference_file, process_metadata, add_region, reorganize_metadata_columns, load_clinical_files, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_lists, merge_wastewater_bams, align_clinical_reads, varmint
except:
    from cli_funcs import process_reference_file, process_metadata, add_region, reorganize_metadata_columns, load_clinical_files, create_monthly_accession_lists, split_clinical_fasta_by_month, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_lists, merge_wastewater_bams, align_clinical_reads, varmint

# entry point function for the CLI
def flu_cli():
    # start timer
    cli_start_time = time.perf_counter()
    # set default file path for output to the current working directory of the user
    output_path_default = os.getcwd()

    # create parser and describe function of script
    parser = argparse.ArgumentParser(description="Process and align wastewater and clinical virome sequencing data for mutation analysis.")

    ############################## set arguments ##############################
    ## required arguments
    # argument for IAV subtype
    parser.add_argument("-s", "--subtype", type=str, required=True, choices=["H1N1", "H3N2", "H5N1"], help="Influenza subtype to process, options are H1N1, H3N2, and H5N1")

    # argument for file path to folder containing wastewater metadata files
    parser.add_argument("-m", "--wastewater_metadata", type=str, required=True, help="Path to folder containing wastewater metadata files")

    # mutually exclusive group for wastewater reads type
    reads_group = parser.add_mutually_exclusive_group(required=True)
    reads_group.add_argument("-pr", "--paired_reads", type=str, help="Path to folders containing paired FASTQ wastewater reads (R1/R2)")
    reads_group.add_argument("-sr", "--single_reads", type=str, help="Path to folder containing single FASTA wastewater reads")

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

    # argument to view current version
    parser.add_argument("-v", "--version", action='store_true', help="View current version")

    # argument to keep all output
    parser.add_argument("-a", "--all", action='store_true', help="Keep all intermediate alignment files created")

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
    
    ############################## process reference and metadata files ##############################
    # optionally give current version
    if args.version:
        logger.info(f"Current version is {version('flu_CLI')}")

    # create directory for unzipped reference files
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    # crm: being too specific about spacing, I know my character flaws
    print("")

    ## process reference files
    section_start = time.perf_counter()
    process_reference_file(args.reference_files, dirs["reference_dir"])
    logger.info(f"Reference processing: {time.perf_counter() - section_start:.2f}s")

    # process wastewater metadata
    section_start = time.perf_counter()
    metadata = process_metadata(args.wastewater_metadata)
    if metadata.empty:
        logger.error("No metadata files found in the specified directory.")
        sys.exit(1)

    # add region column to the metadata if user didn't specify month only
    if not args.month_only:
        # add region column to metadata
        metadata = add_region(metadata)
        no_region = False
    else:
        no_region = True
    
    # reorganize metadata columns with consideration of region potentially not being included
    metadata = reorganize_metadata_columns(metadata, no_region=no_region)
    
    # optionally give time range of the wastewater samples
    if args.time_range:
        earliest, latest = metadata["Date"].min(), metadata["Date"].max()
        logger.info(f"Date range covered by wastewater data: {earliest.strftime('%Y-%m-%d')} to {latest.strftime('%Y-%m-%d')}")

    # create directory for processed and merged metadata
    dirs["metadata_dir"] = os.path.join(dirs["output"], "metadata_files")
    os.makedirs(dirs["metadata_dir"], exist_ok=True)

    # export processed wastewater metadata
    logger.info(f"Exporting the processed metadata to {dirs['metadata_dir']}\n")
    metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_wastewater_combined.csv"), index=False)

    # load in clinical metadata and fasta
    if include_clinical:
        clinical_metadata, clinical_fasta = load_clinical_files(args.clinical_files)

        # export processed clinical metadata
        clinical_metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_clinical_{args.subtype}.csv"), index=False)
    logger.info(f"Metadata processing: {time.perf_counter() - section_start:.2f}s")

    ############################## wastewater ##############################
    # find wastewater reads from pools
    section_start = time.perf_counter()
    # determine which reads type was provided
    if args.single_reads:
        wastewater_reads = find_wastewater_reads(args.single_reads, args.subtype, single_reads=True)
    else:
        wastewater_reads = find_wastewater_reads(args.paired_reads, args.subtype, single_reads=False)
    logger.info(f"Finding wastewater reads: {time.perf_counter() - section_start:.2f}s")
    
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
    section_start = time.perf_counter()
    align_wastewater_reads(wastewater_reads, dirs["reference_dir"], dirs["pools"])
    logger.info(f"\nWastewater alignment: {time.perf_counter() - section_start:.2f}s")

    # create directory for wastewater lists
    dirs["wastewater_lists_dir"] = os.path.join(dirs["wastewater_dir"], "lists")
    os.makedirs(dirs["wastewater_lists_dir"], exist_ok=True)

    # create subfolders for wastewater lists
    if include_region:
        dirs["wastewater_list_month"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month")
        os.makedirs(dirs["wastewater_list_month"], exist_ok=True)
        dirs["wastewater_list_region"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month_region")
        os.makedirs(dirs["wastewater_list_region"], exist_ok=True)
    else:
        dirs["wastewater_list_month"] = dirs["wastewater_lists_dir"]

    # create list of monthly accessions for downstream merge key creation
    if include_region:
        logger.info("\nMerging wastewater alignment files by month and public health region")
    else:
        logger.info("\nMerging wastewater alignment files by month")

    section_start = time.perf_counter()
    create_wastewater_bam_lists(metadata, dirs["pools"], dirs["wastewater_list_month"], dirs.get("wastewater_list_region"), include_region)
    logger.info(f"Creating BAM lists: {time.perf_counter() - section_start:.2f}s")

    # wastewater merged bams
    dirs["merged_bams"] = os.path.join(dirs["wastewater_dir"], "merged_bams")
    os.makedirs(dirs["merged_bams"], exist_ok=True)

    # create subfolder: wastewater bams merged by month
    dirs["merged_bams_month"] = os.path.join(dirs["merged_bams"], "merged_bams_month")        
    os.makedirs(dirs["merged_bams_month"], exist_ok=True)

    # merge wastewater bams by month
    section_start = time.perf_counter()
    merge_wastewater_bams(dirs["wastewater_list_month"], dirs["merged_bams_month"])
    logger.info(f"Merging BAMs by month: {time.perf_counter() - section_start:.2f}s")

    # create subfolder: wastewater bams merged by month and region
    if include_region:
        dirs["merged_bams_month_region"] = os.path.join(dirs["merged_bams"], "merged_bams_month_region")
        os.makedirs(dirs["merged_bams_month_region"], exist_ok=True)
        
        # merge wastewater bams by month and region
        section_start = time.perf_counter()
        merge_wastewater_bams(dirs["wastewater_list_region"], dirs["merged_bams_month_region"])
        logger.info(f"Merging BAMs by month+region: {time.perf_counter() - section_start:.2f}s")
    
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
            
    # varmint for wastewater
    logger.info("\nAnnotating coding effects of mutations with varmint")
    section_start = time.perf_counter()

    if include_region:
        varmint(dirs["merged_bams_month"], dirs["reference_dir"], dirs["tsv_month"])
        varmint(dirs["merged_bams_month_region"], dirs["reference_dir"], dirs["tsv_month_region"])
    else:
        varmint(dirs["merged_bams_month"], dirs["reference_dir"], dirs["tsv_output"])
    
    logger.info(f"Varmint (wastewater): {time.perf_counter() - section_start:.2f}s")

    ############################## clinical ##############################
    if include_clinical:
        # create folder for clinical output
        dirs["clinical"] = os.path.join(dirs["alignment_dir"], "clinical")
        os.makedirs(dirs["clinical"], exist_ok=True)

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
        section_start = time.perf_counter()
        split_clinical_fasta_by_month(clinical_fasta, dirs["clinical_lists_month"], dirs["clinical_fasta_month"])
        logger.info(f"Splitting clinical FASTA: {time.perf_counter() - section_start:.2f}s")

        # create folder for the clinical bam files that were merged by month
        dirs["clinical_bam_month"] = os.path.join(dirs["clinical"], "clinical_bam_month")
        os.makedirs(dirs["clinical_bam_month"], exist_ok=True)

        # align clinical reads to reference
        logger.info("\nAligning clinical reads to the reference genome")
        section_start = time.perf_counter()
        align_clinical_reads(dirs["clinical_fasta_month"], dirs["clinical_bam_month"], dirs["reference_dir"])
        logger.info(f"Clinical alignment: {time.perf_counter() - section_start:.2f}s")

        # varmint for clinical
        logger.info("\nAnnotating coding effects of mutations with varmint")
        section_start = time.perf_counter()
        varmint(dirs["clinical_bam_month"], dirs["reference_dir"], dirs["tsv_clinical"])
        logger.info(f"Varmint (clinical): {time.perf_counter() - section_start:.2f}s")

        # crm: maybe doing too much here
        # delete alignment files if not requested
        if not args.all:
            shutil.rmtree(dirs["alignment_dir"])
            logger.info(f"\nRemoved intermediate alignment files")

    # crm: check iav_serotype and make sure you're ending this correctly
    # print run time
    cli_end_time = time.perf_counter()
    time_taken = round((cli_end_time - cli_start_time), 2) 
    # crm: can make this prettier later
    logger.info(f"Run time: {timedelta(seconds=time_taken)}")