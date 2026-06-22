#!/usr/bin/env python

###################### SETUP ######################
# load modules
import argparse
import sys, os
import logging
import pysam
import shutil
import time                         
from datetime import timedelta      

# CPU detection for fast mode
cpu_count = os.cpu_count() or 4

# load in functions from mutatea_covid_funcs
try:
    from .mutatea_covid_funcs import (process_reference_files, process_metadata, add_region, create_wastewater_bam_groups, merge_wastewater_bams, 
    find_bam_files_from_covid_samples, run_stats, varmint, print_mutatea_banner)
except:
    from mutatea_covid_funcs import (process_reference_files, process_metadata, add_region, create_wastewater_bam_groups, merge_wastewater_bams, 
    find_bam_files_from_covid_samples, run_stats, varmint, print_mutatea_banner)

# entry point function for the CLI
def mutatea_covid():
    # print ASCII art over time
    for line in print_mutatea_banner().splitlines():
        print(line)
        time.sleep(0.04)

    # start timer
    cli_start_time = time.perf_counter()
    # set default file path for output to the current working directory of the user
    output_path_default = os.getcwd()

    # create parser and describe function of script
    parser = argparse.ArgumentParser(description="Process wastewater BAM files for mutation analysis - COVID positive samples version.")

    ############################## set arguments ##############################
    ## required arguments
    # argument for file path to folder containing wastewater metadata files
    parser.add_argument("-m", "--wastewater_metadata", type=str, required=True, help="Path to folder containing wastewater metadata files")

    # argument for file path to covid positive samples file
    parser.add_argument("-cps", "--covid_positive_samples", type=str, required=True, help="Path to file containing COVID positive samples for targeted analysis")

    # argument for wastewater BAM files
    parser.add_argument("-b", "--bam_files", type=str, required=True, help="Path to directory containing wastewater BAM files")

    # argument for file path to folder containing reference files
    parser.add_argument("-ref", "--reference_files", type=str, required=True, help="Path to folder containing the reference fasta(.gz) and gff(.gz) files")

    ## optional arguments
    
    # argument for nondefault output directory
    parser.add_argument("-o", "--output_dir", type=str, default=output_path_default, help="Path to chosen output directory")

    # argument to only split wastewater metadata by time
    parser.add_argument("-ty", "--time_only", action='store_true', help="Override default split of time and region, will only split the wastewater data by time")

    # argument to view time range covered by wastewater metadata
    parser.add_argument("-tr", "--timerange", action='store_true', help="View time range covered by wastewater sample collection")

    # argument to view current version
    parser.add_argument("-v", "--version", action='store_true', help="View current version")

    # argument to keep all output
    parser.add_argument("-a", "--all", action='store_true', help="Keep all intermediate files")

    # argument to override default number for parallel workers
    parser.add_argument("-f", "--fast", action='store_true', help="Override default number of parallel workers to run with all available cpus")

    # argument to input personal dictionary (default is mapping Texas city to Texas public health region)
    parser.add_argument("-d", "--dictionary", type=str, help="Path to JSON file containing city-to-region mapping")

    # argument to save detailed loggger file
    parser.add_argument("-l", "--logger", action='store_true', help="Export a detailed logger file")

    # argument to save statistics of the groupings
    parser.add_argument("-s", "--statistics", action='store_true', help="Export a file detailing the genome depth and coverage for each grouping") 
    
    # argument to change the unit of time by which the samples are organized (default is month)
    parser.add_argument("-g", "--grouping", choices=["month", "week", "day", "year"], help="Group samples by year, month, week, or day (default is month)")

    # argument to overwrite default mapq filter (default is 0)
    parser.add_argument("-q", "--mapq", type=int, default=0, help="Minimum mapping quality (MAPQ) threshold for filtering reads (default is 0)")
    
    # parse arguments
    args = parser.parse_args()

    # ensure grouping is lowercase and without spaces
    if args.grouping:
        args.grouping = args.grouping.lower().replace(" ", "")
    grouping = args.grouping or "month"

    ## boolean checks
    # check if region is included
    include_region = not args.time_only

    # validate covid positive samples file
    if not os.path.exists(args.covid_positive_samples):
        sys.exit(f"Error: COVID positive samples file not found: {args.covid_positive_samples}")

    # initialize directories dictionary
    dirs = {}
    
    # create main output directory
    dirs["output"] = os.path.join(args.output_dir, "covid_align")
    os.makedirs(dirs["output"], exist_ok=True)

    # define logger
    logger = logging.getLogger("mutatea_covid_logger")
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)

    logger.addHandler(stream_handler)

    # save detailed log file if user requested
    if args.logger:
        log_file = os.path.join(dirs["output"], "covid_mutatea.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(file_handler)
    
    ############################## process reference and metadata files ##############################
    # optionally give current version
    if args.version:
        logger.info(f"Current version is COVID-positive targeted analysis")

    # create directory for unzipped reference files
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    print("")

    ## process reference files
    section_start = time.perf_counter()
    try:
        fna_path, gff_path = process_reference_files(args.reference_files, dirs["reference_dir"])
    except Exception as e:
        sys.exit(f"Error processing the reference files: {e}")  
    logger.info(f"Reference processing: {time.perf_counter() - section_start:.2f}s")

    # process wastewater metadata
    section_start = time.perf_counter()
    try:
        metadata = process_metadata(args.wastewater_metadata, grouping, logger=logger)
    except Exception as e:
        sys.exit(f"Error processing metadata: {e}")
    
    if metadata.empty:
        logger.error("No metadata files found in the specified directory.")
        sys.exit(1)

    # add region column to the metadata if user didn't specify time_only grouping
    if not args.time_only:
        try:
            # add region column to metadata
            metadata = add_region(metadata, args.dictionary)
        except Exception as e:
            sys.exit(f"Error adding region to the metadata: {e}") 

    # optionally give time range of the wastewater samples
    if args.timerange:
        earliest, latest = metadata["Date"].min(), metadata["Date"].max()
        logger.info(f"Date range covered by wastewater data: {earliest.strftime('%Y-%m-%d')} to {latest.strftime('%Y-%m-%d')}")

    # create directory for processed and merged metadata
    dirs["metadata_dir"] = os.path.join(dirs["output"], "metadata_files")
    os.makedirs(dirs["metadata_dir"], exist_ok=True)

    # export processed wastewater metadata
    logger.info(f"Exporting the processed metadata to {dirs['metadata_dir']}")
    try:
        metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_wastewater_combined.csv"), index=False)
    except Exception as e:
        sys.exit(f"Error exporting the processed metadata: {e}") 

    logger.info(f"Metadata processing: {time.perf_counter() - section_start:.2f}s")

    ############################## wastewater ##############################
    # find wastewater reads from covid positive samples
    section_start = time.perf_counter()
    print("")
    logger.info(f"Finding wastewater reads from COVID positive samples")
    
    # find BAM files for COVID positive samples
    try:
        wastewater_bams = find_bam_files_from_covid_samples(args.covid_positive_samples, args.bam_files)
    except Exception as e:
        sys.exit(f"Error finding the wastewater BAM files from COVID samples: {e}") 
    logger.info(f"Finding reads (wastewater): {time.perf_counter() - section_start:.2f}s")
    
    # filter BAM files by MAPQ if specified
    bam_files = wastewater_bams
    if args.mapq > 0:
        logger.info(f"Filtering BAM files by MAPQ >= {args.mapq}")
        filtered_bams = []
        for bam_file in bam_files:
            filtered_bam = bam_file.replace(".bam", f".mapq{args.mapq}.bam")
            cmd_filter = f"samtools view -b -q {args.mapq} {bam_file} > {filtered_bam}"
            try:
                subprocess.run(cmd_filter, shell=True, check=True, executable="/bin/bash")
                filtered_bams.append(filtered_bam)
            except Exception as e:
                logger.warning(f"Could not filter {bam_file}: {e}")
                filtered_bams.append(bam_file)
        bam_files = filtered_bams
    
    logger.info(f"Processing {len(bam_files)} BAM files from COVID-positive pools")

    # create directory for wastewater processing
    dirs["wastewater_dir"] = os.path.join(dirs["output"], "wastewater")
    os.makedirs(dirs["wastewater_dir"], exist_ok=True)

    # create directory for wastewater lists
    dirs["wastewater_lists_dir"] = os.path.join(dirs["wastewater_dir"], "lists")
    os.makedirs(dirs["wastewater_lists_dir"], exist_ok=True)

    # create subfolders for wastewater lists
    if include_region:
        dirs[f"wastewater_list_{grouping}"] = os.path.join(dirs["wastewater_lists_dir"], f"lists_{grouping}")
        os.makedirs(dirs[f"wastewater_list_{grouping}"], exist_ok=True)
        dirs["wastewater_list_region"] = os.path.join(dirs["wastewater_lists_dir"], f"lists_{grouping}_region")
        os.makedirs(dirs["wastewater_list_region"], exist_ok=True)
    else:
        dirs[f"wastewater_list_{grouping}"] = dirs["wastewater_lists_dir"]

    # create list of accessions by time grouping for downstream merge key creation
    print("")
    if include_region:
        logger.info(f"Merging wastewater alignment files by {grouping} and region")
    else:
        logger.info(f"Merging wastewater alignment files by {grouping}")

    section_start = time.perf_counter()
    try:
        if include_region:
            month_list_dir, region_list_dir = create_wastewater_bam_groups(bam_files, metadata, dirs[f"wastewater_list_{grouping}"], dirs.get("wastewater_list_region"), include_region, grouping)
        else:
            month_list_dir = create_wastewater_bam_groups(bam_files, metadata, dirs[f"wastewater_list_{grouping}"], dirs.get("wastewater_list_region"), include_region, grouping)
    except Exception as e:
        sys.exit(f"Error creating the lists for merging wastewater alignment files: {e}") 
    logger.info(f"Creating BAM lists: {time.perf_counter() - section_start:.2f}s")

    # wastewater merged bams
    merged_dir_name = "bams_filtered_merged" if args.mapq else "bams_merged"
    dirs["bams_merged"] = os.path.join(dirs["wastewater_dir"], merged_dir_name)
    os.makedirs(dirs["bams_merged"], exist_ok=True)

    # create subfolder: wastewater bams merged by chosen time grouping
    dirs[f"bams_{grouping}"] = os.path.join(dirs["bams_merged"], f"bams_{grouping}")        
    os.makedirs(dirs[f"bams_{grouping}"], exist_ok=True)

    # merge wastewater bams by chosen time grouping
    section_start = time.perf_counter()
    try:
        merged_bams_time = merge_wastewater_bams(month_list_dir, dirs[f"bams_{grouping}"], min_mapq=args.mapq)
    except Exception as e:
        sys.exit(f"Error merging wastewater alignment files by {grouping}: {e}") 
    logger.info(f"Merging BAMs by {grouping}: {time.perf_counter() - section_start:.2f}s")

    # create subfolder: wastewater bams merged by chosen time grouping and region
    if include_region:
        dirs[f"bams_{grouping}_region"] = os.path.join(dirs["bams_merged"], f"bams_{grouping}_region")
        os.makedirs(dirs[f"bams_{grouping}_region"], exist_ok=True)
        
        # merge wastewater bams by chosen time grouping and region
        section_start = time.perf_counter()
        try:
            merged_bams_time_region = merge_wastewater_bams(region_list_dir, dirs[f"bams_{grouping}_region"], min_mapq=args.mapq)
        except Exception as e:
            sys.exit(f"Error creating the lists for merging wastewater alignment files by chose time grouping and region: {e}") 
        logger.info(f"Merging BAMs by {grouping} and region: {time.perf_counter() - section_start:.2f}s")
    
    # get genome coverage if statistics included
    if args.statistics:
        print("")
        logger.info("Getting coverage statistics of wastewater BAMs with samtools")
        section_start = time.perf_counter()

        # create stats folder to later catch tsv files
        dirs["statistics"] = os.path.join(dirs["output"], "statistics")
        os.makedirs(dirs["statistics"], exist_ok=True)

        if include_region:
            dirs[f"statistics_{grouping}"] = os.path.join(dirs["statistics"], f"statistics_{grouping}")
            os.makedirs(dirs[f"statistics_{grouping}"], exist_ok=True)
            dirs[f"statistics_{grouping}_region"] = os.path.join(dirs["statistics"], f"statistics_{grouping}_region")
            os.makedirs(dirs[f"statistics_{grouping}_region"], exist_ok=True)
            
            try:
                statistics = run_stats(merged_bams_time, dirs[f"statistics_{grouping}"], logger=logger)
                statistics = run_stats(merged_bams_time_region, dirs[f"statistics_{grouping}_region"], logger=logger)
            except Exception as e:
                sys.exit(f"Error getting coverage statistics of wastewater BAMs with samtools: {e}")
        else:
            dirs[f"statistics_{grouping}"] = dirs["statistics"]
            
            try:
                statistics = run_stats(merged_bams_time, dirs["statistics"], logger=logger)
            except Exception as e:
                sys.exit(f"Error getting coverage statistics of wastewater BAMs with samtools: {e}")
        logger.info(f"Running statistics on merged BAMs (wastewater): {time.perf_counter() - section_start:.2f}s")

    # create tsv_output folder to later catch tsv files
    dirs["tsv_output"] = os.path.join(dirs["output"], "tsv_output")
    os.makedirs(dirs["tsv_output"], exist_ok=True)

    # create subfolders for TSV output
    if include_region:
        dirs[f"tsv_{grouping}"] = os.path.join(dirs["tsv_output"], f"{grouping}")
        os.makedirs(dirs[f"tsv_{grouping}"], exist_ok=True)
        dirs[f"tsv_{grouping}_region"] = os.path.join(dirs["tsv_output"], f"{grouping}_region")
        os.makedirs(dirs[f"tsv_{grouping}_region"], exist_ok=True)
    else:
        dirs[f"tsv_{grouping}"] = dirs["tsv_output"]
    
    # Run varmint
    print("")
    logger.info("Annotating coding effects of mutations with varmint (wastewater samples)")
    section_start = time.perf_counter()
    
    if include_region:
        try:
            print(f"Running varmint (time) using {cpu_count if args.fast else 4} parallel workers")
            varmint(merged_bams_time, fna_path, gff_path, dirs[f"tsv_{grouping}"], workers=cpu_count if args.fast else 4)
            print(f"Running varmint (time+region) using {cpu_count if args.fast else 4} parallel workers")
            varmint(merged_bams_time_region, fna_path, gff_path, dirs[f"tsv_{grouping}_region"], workers=cpu_count if args.fast else 4)
        except Exception as e:
            sys.exit(f"Error running varmint: {e}")
    else:
        try:
            print(f"Running varmint (time) using {cpu_count if args.fast else 4} parallel workers")
            varmint(merged_bams_time, fna_path, gff_path, dirs["tsv_output"], workers=cpu_count if args.fast else 4)
        except Exception as e:
            sys.exit(f"Error running varmint: {e}")
    
    logger.info(f"Varmint (wastewater): {time.perf_counter() - section_start:.2f}s")


    # Note: Clinical processing is not supported in this BAM-only version
    # Clinical data processing requires FASTA files and alignment, which are not handled here

    # print run time
    cli_end_time = time.perf_counter()
    time_taken = round((cli_end_time - cli_start_time), 2) 
    print("")
    logger.info(f"Run time of mutatea_covid: {timedelta(seconds=time_taken)}")

if __name__ == "__main__":
    mutatea_covid()
