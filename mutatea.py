#!/usr/bin/env python

###################### SETUP ######################
# load modules
import argparse
import sys, os
import logging
import shutil
import time                         # used for granular timing of each function called in
from datetime import timedelta      # used to determine total runtime of mutatea run

# CPU detection for fast mode
cpu_count = os.cpu_count() or 4

# load in functions from mutatea.funcs
try:
    from .mutatea_funcs import process_reference_files, process_metadata, add_region, load_clinical_files, create_grouped_accession_lists, split_clinical_fasta_by_time, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_groups, merge_wastewater_bams, align_clinical_reads, run_stats, varmint, alignment_quality_filter
except:
    from mutatea_funcs import process_reference_files, process_metadata, add_region, load_clinical_files, create_grouped_accession_lists, split_clinical_fasta_by_time, find_wastewater_reads, align_wastewater_reads, create_wastewater_bam_groups, merge_wastewater_bams, align_clinical_reads, run_stats, varmint, alignment_quality_filter

# entry point function for the CLI
def mutatea():
    # start timer
    cli_start_time = time.perf_counter()
    # set default file path for output to the current working directory of the user
    output_path_default = os.getcwd()

    # create parser and describe function of script
    parser = argparse.ArgumentParser(description="Process and align wastewater and clinical pathogen sequencing data for mutation analysis.")

    ############################## set arguments ##############################
    ## required arguments
    # argument for pathogen
    parser.add_argument("-p", "--pathogen", type=str, required=True, help="Pathogen being processed, name should match the naming of the reads")

    # argument for file path to folder containing wastewater metadata files
    parser.add_argument("-m", "--wastewater_metadata", type=str, required=True, help="Path to folder containing wastewater metadata files")

    # mutually exclusive group for wastewater reads type
    reads_group = parser.add_mutually_exclusive_group(required=True)
    reads_group.add_argument("-pr", "--paired_reads", type=str, help="Path to folders containing paired FASTQ wastewater reads (R1/R2)")
    reads_group.add_argument("-sr", "--single_reads", type=str, help="Path to folder containing single FASTQ wastewater reads")

    # argument for file path to folder containing reference files
    parser.add_argument("-ref", "--reference_files", type=str, required=True, help="Path to folder containing the reference fasta(.gz) and gff(.gz) files")

    ## optional arguments
    # argument for file path to folder containing clinical files
    parser.add_argument("-c", "--clinical_files", type=str, help="Path to folder containing the clinical fasta and csv files")

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

    # argument to input personal dictionary (default is mapping Texas city to public health region)
    parser.add_argument("-d", "--dictionary", type=str, help="Path to JSON file containing city-to-region mapping")

    # argument to save detailed loggger file
    parser.add_argument("-l", "--logger", action='store_true', help="Export a detailed logger file")

    # crm: wants to add in the statistics argument
    # argument to save statistics of the groupings
    parser.add_argument("-s", "--statistics", action='store_true', help="Export a file detailing the genome depth and coverage for each grouping") 
    
    # argument to change the unit of time by which the samples are organized (default is month)
    parser.add_argument("-g", "--grouping", choices=["month", "week", "day", "year"], help="Group samples by year, month, week, or day (default is month)")

    # argument to overwrite minimap2 preset for the wastewater alignment (default is -ax sr)
    parser.add_argument("-mw", "--minimap_wastewater", type=str, default="sr", help="Override minimap2 preset for wastewater alignment (default is sr)")

    # argument to overwrite minimap2 preset for the clinical alignment (default is -ax asm10)
    parser.add_argument("-mc", "--minimap_clinical", type=str, default="asm10", help="Override minimap2 preset for clinical alignment (default is asm10)")

    # parse arguments
    args = parser.parse_args()

    # ensure grouping is lowercase and without spaces
    if args.grouping:
        args.grouping = args.grouping.lower().replace(" ", "")
    grouping = args.grouping or "month"

    ## boolean checks
    # check if clinical files are included
    include_clinical = bool(args.clinical_files)

    # check if region is included
    include_region = not args.time_only

    # initialize directories dictionary
    dirs = {}
    
    # create main output directory with pathogen-specific subfolder
    dirs["output"] = os.path.join(args.output_dir, f"{args.pathogen}_align")
    os.makedirs(dirs["output"], exist_ok=True)

    # crm: maybe add more here
    # define logger
    logger = logging.getLogger("mutatea_logger")
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)

    logger.addHandler(stream_handler)

    # save detailed log file if user requested
    if args.logger:
        log_file = os.path.join(dirs["output"], f"{args.pathogen}_mutatea.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(file_handler)
    
    ############################## process reference and metadata files ##############################
    # optionally give current version
    if args.version:
        logger.info(f"Current version is {version('mutatea')}")

    # create directory for unzipped reference files
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    # crm: being too specific about spacing, I know my character flaws
    logger.info("")

    ## process reference files
    section_start = time.perf_counter()
    try:
        fna_path, gff_path = process_reference_files(args.reference_files, dirs["reference_dir"])
    except Exception as e:
        return f"Error processing the reference files: {e}"  
    logger.info(f"Reference processing: {time.perf_counter() - section_start:.2f}s\n")

    # process wastewater metadata
    section_start = time.perf_counter()
    try:
        metadata = process_metadata(args.wastewater_metadata, getattr(args, 'grouping', 'month'))
    except Exception as e:
        return f"Error processing metadata: {e}"
    
    if metadata.empty:
        logger.error("No metadata files found in the specified directory.")
        sys.exit(1)

    # add region column to the metadata if user didn't specify time_only grouping
    if not args.time_only:
        try:
            # add region column to metadata
            metadata = add_region(metadata, args.dictionary)
        except Exception as e:
            return f"Error adding region to the metadata: {e}" 

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
        return f"Error exporting the processed metadata: {e}" 

    # load in clinical metadata and fasta
    if include_clinical:
        try:
            clinical_metadata, clinical_fasta = load_clinical_files(args.clinical_files, args.grouping)

            # export processed clinical metadata
            clinical_metadata.to_csv(os.path.join(dirs["metadata_dir"], f"metadata_clinical_{args.pathogen}.csv"), index=False)
        except Exception as e:
            return f"Error loading in the clinical files: {e}" 
    logger.info(f"Metadata processing: {time.perf_counter() - section_start:.2f}s\n")

    ############################## wastewater ##############################
    # find wastewater reads from pools
    section_start = time.perf_counter()
    logger.info(f"Finding wastewater reads from pools")

    # determine which read type was provided
    if args.single_reads:
        try:
            wastewater_reads = find_wastewater_reads(args.single_reads, args.pathogen, single_reads=True)
        except Exception as e:
            return f"Error finding the wastewater reads: {e}" 
    else:
        try:
            wastewater_reads = find_wastewater_reads(args.paired_reads, args.pathogen, single_reads=False)
        except Exception as e:
            return f"Error finding the wastewater reads: {e}" 
    logger.info(f"Finding reads (wastewater): {time.perf_counter() - section_start:.2f}s")
    
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
    logger.info("\nAligning wastewater reads to given reference genome")
    section_start = time.perf_counter()
    try:
        bam_files = align_wastewater_reads(wastewater_reads, fna_path, dirs["pools"], pathogen=args.pathogen, minimap_preset=args.minimap_wastewater, workers=cpu_count if args.fast else 4)
    except Exception as e:
        return f"Error aligning the wastewater reads: {e}" 
    logger.info(f"Aligning reads to reference genome (wastewater): {time.perf_counter() - section_start:.2f}s")

    # create directory for mapq filtered reads
    dirs["wastewater_filtered"] = os.path.join(dirs["wastewater_dir"], "filtered_merged_bams")
    os.makedirs(dirs["wastewater_filtered"], exist_ok=True)

    # filter bams for mapping quality
    logger.info(f"\nFiltering wastewater reads for mapping quality")
    section_start = time.perf_counter()
    try:
        bam_files = alignment_quality_filter(bam_files, dirs["wastewater_filtered"])
    except Exception as e:
        return f"Error filtering wastewater reads for mapping quality: {e}" 
    logger.info(f"Filtering reads for mapping quality (wastewater): {time.perf_counter() - section_start:.2f}s")


    # crm: need to add in a print line explaining how many reads were removed and from where (also for clinical)

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
    if include_region:
        logger.info(f"\nMerging wastewater alignment files by {grouping} and public health region")
    else:
        logger.info(f"\nMerging wastewater alignment files by {grouping}")

    section_start = time.perf_counter()
    try:
        if include_region:
            month_list_dir, region_list_dir = create_wastewater_bam_groups(bam_files, metadata, dirs[f"wastewater_list_{grouping}"], dirs.get("wastewater_list_region"), include_region)
        else:
            month_list_dir = create_wastewater_bam_groups(bam_files, metadata, dirs[f"wastewater_list_{grouping}"], dirs.get("wastewater_list_region"), include_region)
    except Exception as e:
        return f"Error creating the lists for merging wastewater alignment files: {e}" 
    logger.info(f"Creating BAM lists: {time.perf_counter() - section_start:.2f}s")

    # wastewater merged bams
    dirs["merged_bams"] = os.path.join(dirs["wastewater_dir"], "merged_bams")
    os.makedirs(dirs["merged_bams"], exist_ok=True)

    # create subfolder: wastewater bams merged by chosen time grouping
    dirs[f"merged_bams_{grouping}"] = os.path.join(dirs["merged_bams"], f"merged_bams_{grouping}")        
    os.makedirs(dirs[f"merged_bams_{grouping}"], exist_ok=True)

    # merge wastewater bams by chosen time grouping
    section_start = time.perf_counter()
    try:
        merged_bams_time = merge_wastewater_bams(month_list_dir, dirs[f"merged_bams_{grouping}"])
    except Exception as e:
        return f"Error merging wastewater alignment files by {grouping}: {e}" 
    logger.info(f"Merging BAMs by {grouping}: {time.perf_counter() - section_start:.2f}s")

    # create subfolder: wastewater bams merged by chosen time grouping and region
    if include_region:
        dirs[f"merged_bams_{grouping}_region"] = os.path.join(dirs["merged_bams"], f"merged_bams_{grouping}_region")
        os.makedirs(dirs[f"merged_bams_{grouping}_region"], exist_ok=True)
        
        # merge wastewater bams by chosen time grouping and region
        section_start = time.perf_counter()
        try:
            merged_bams_time_region = merge_wastewater_bams(region_list_dir, dirs[f"merged_bams_{grouping}_region"])
        except Exception as e:
            return f"Error creating the lists for merging wastewater alignment files by chose time grouping and region: {e}" 
        logger.info(f"Merging BAMs by {grouping} and region: {time.perf_counter() - section_start:.2f}s")
    
    # get genome coverage if statistics included
    if args.statistics:
        logger.info("\nGetting coverage statistics of wastewater BAMs with samtools")
        section_start = time.perf_counter()

        # create stats folder to later catch tsv files
        dirs["statistics"] = os.path.join(dirs["output"], "statistics")
        os.makedirs(dirs["statistics"], exist_ok=True)

        # create subfolders if clinical included
        if include_clinical:
            # split output statistics files by source
            dirs["stats_wastewater"] = os.path.join(dirs["statistics"], "wastewater")
            os.makedirs(dirs["stats_wastewater"], exist_ok=True)
            dirs["stats_clinical"] = os.path.join(dirs["statistics"], "clinical")
            os.makedirs(dirs["stats_clinical"], exist_ok=True)
            
            # split wastewater output by grouping method
            if include_region:
                # create subfolders
                dirs[f"statistics_{grouping}"] = os.path.join(dirs["stats_wastewater"], f"statistics_{grouping}")
                os.makedirs(dirs[f"statistics_{grouping}"], exist_ok=True)
                dirs[f"statistics_{grouping}_region"] = os.path.join(dirs["stats_wastewater"], f"statistics_{grouping}_region")
                os.makedirs(dirs[f"statistics_{grouping}_region"], exist_ok=True)
                
                try:
                    statistics = run_stats(merged_bams_time, dirs[f"statistics_{grouping}"])
                    statistics = run_stats(merged_bams_time_region, dirs[f"statistics_{grouping}_region"])
                except Exception as e:
                    return f"Error getting coverage statistics of wastewater BAMs with samtools: {e}" 
            else:
                try:
                    statistics = run_stats(merged_bams_time, dirs["stats_wastewater"])
                except Exception as e:
                    return f"Error getting coverage statistics of wastewater BAMs with samtools: {e}"

        else:
            if include_region:
                dirs[f"statistics_{grouping}"] = os.path.join(dirs["statistics"], f"statistics_{grouping}")
                os.makedirs(dirs[f"statistics_{grouping}"], exist_ok=True)
                dirs[f"statistics_{grouping}_region"] = os.path.join(dirs["statistics"], f"statistics_{grouping}_region")
                os.makedirs(dirs[f"statistics_{grouping}_region"], exist_ok=True)
                
                try:
                    statistics = run_stats(merged_bams_time, dirs[f"statistics_{grouping}"])
                    statistics = run_stats(merged_bams_time_region, dirs[f"statistics_{grouping}_region"])
                except Exception as e:
                    return f"Error getting coverage statistics of wastewater BAMs with samtools: {e}" 
            else:
                try:
                    statistics = run_stats(merged_bams_time, dirs["statistics"])
                except Exception as e:
                    return f"Error getting coverage statistics of wastewater BAMs with samtools: {e}"
        logger.info(f"Running statistics on merged BAMs (wastewater): {time.perf_counter() - section_start:.2f}s")

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
            dirs[f"tsv_{grouping}"] = os.path.join(dirs["tsv_wastewater"], f"{grouping}")
            os.makedirs(dirs[f"tsv_{grouping}"], exist_ok=True)
            dirs[f"tsv_{grouping}_region"] = os.path.join(dirs["tsv_wastewater"], f"{grouping}_region")
            os.makedirs(dirs[f"tsv_{grouping}_region"], exist_ok=True)

    # or just split clinical output by grouping method
    else:
        dirs[f"tsv_{grouping}"] = os.path.join(dirs["tsv_output"], f"{grouping}")
        os.makedirs(dirs[f"tsv_{grouping}"], exist_ok=True)

        dirs[f"tsv_{grouping}_region"] = os.path.join(dirs["tsv_output"], f"{grouping}_region")
        os.makedirs(dirs[f"tsv_{grouping}_region"], exist_ok=True)
    
    # Run varmint
    logger.info("\nAnnotating coding effects of mutations with varmint (wastewater samples)")
    section_start = time.perf_counter()
    
    if include_region:
        try:
            varmint(merged_bams_time, fna_path, gff_path, dirs[f"tsv_{grouping}"], workers=cpu_count if args.fast else 4)
            varmint(merged_bams_time_region, fna_path, gff_path, dirs[f"tsv_{grouping}_region"], workers=cpu_count if args.fast else 4)
        except Exception as e:
            return f"Error running varmint: {e}"
    else:
        try:
            varmint(merged_bams_time, fna_path, gff_path, dirs["tsv_output"], workers=cpu_count if args.fast else 4)
        except Exception as e:
            return f"Error running varmint: {e}"
    
    logger.info(f"Varmint (wastewater): {time.perf_counter() - section_start:.2f}s")

    ############################## clinical ##############################
    if include_clinical:
        # create folder for clinical output
        dirs["clinical"] = os.path.join(dirs["alignment_dir"], "clinical")
        os.makedirs(dirs["clinical"], exist_ok=True)

        # create folder for the lists of accessions by chosen time grouping
        dirs[f"lists_{grouping}"] = os.path.join(dirs["clinical"], f"lists_{grouping}")
        os.makedirs(dirs[f"lists_{grouping}"], exist_ok=True)

        # create monthly lists of accessions
        # crm: I want to delete this logger line, which I could go forward with if I am saving the monthly lists to tempdir
        #logger.info("\nCreating monthly lists of accessions")
        try:
            create_grouped_accession_lists(clinical_metadata, dirs[f"lists_{grouping}"])
        except Exception as e:
            return f"Error creating lists of accessions grouped by time for clinical data: {e}" 

        # create folder for the clinical fastas split by chosen time grouping
        dirs[f"fastas_{grouping}"] = os.path.join(dirs["clinical"], f"fastas_{grouping}")
        os.makedirs(dirs[f"fastas_{grouping}"], exist_ok=True)

        # split clinical fasta by time-grouped lists
        logger.info(f"\nSplitting clinical FASTA by {grouping}")
        section_start = time.perf_counter()
        try:
            split_clinical_fasta_by_time(clinical_fasta, dirs[f"lists_{grouping}"], dirs[f"fastas_{grouping}"])
        except Exception as e:
            return f"Error splitting clinical FASTA by {grouping}: {e}" 
        logger.info(f"Splitting FASTA (clinical): {time.perf_counter() - section_start:.2f}s")

        # create folder for the clinical bam files that were merged by chosen time grouping
        dirs[f"bam_{grouping}"] = os.path.join(dirs["clinical"], f"bam_{grouping}")
        os.makedirs(dirs[f"bam_{grouping}"], exist_ok=True)

        # align clinical reads to reference
        logger.info("\nAligning clinical reads to the reference genome")
        section_start = time.perf_counter()
        try:
            bam_files = align_clinical_reads(dirs[f"fastas_{grouping}"], fna_path, dirs[f"bam_{grouping}"], minimap_preset=args.minimap_clinical, workers=cpu_count if args.fast else 4, grouping=args.grouping)
        except Exception as e:
            return f"Error aligning the clinical reads: {e}"  
        logger.info(f"Aligning reads to reference genome (clinical): {time.perf_counter() - section_start:.2f}s")

        # create folder for the filtered clinical bam files that were merged by chosen time grouping
        dirs[f"bams_{grouping}_filtered"] = os.path.join(dirs[f"bams_{grouping}"], f"bams_{grouping}_filtered")
        os.makedirs(dirs[f"bams_{grouping}_filtered"], exist_ok=True)

        # filter clinical reads for alignment quality
        logger.info("\nFiltering clinical reads for mapping quality")
        section_start = time.perf_counter()
        try:
            bam_files = alignment_quality_filter(bam_files, dirs[f"bams_{grouping}_filtered"])
        except Exception as e:
            return f"Error filtering clinical reads for mapping quality: {e}" 
        logger.info(f"Filtering reads for mapping quality (clinical): {time.perf_counter() - section_start:.2f}s")











        # get genome coverage if statistics included
        if args.statistics:
            # run statistics
            logger.info("\nGetting coverage statistics of clinical BAMs with samtools")
            section_start = time.perf_counter()
            try:
                statistics = run_stats(bam_files, dirs["stats_clinical"])
            except Exception as e:
                return f"Error getting coverage statistics of clinical BAMs with samtools: {e}"
            logger.info(f"Creating coverage statistics from BAMs (clinical): {time.perf_counter() - section_start:.2f}s")

        # varmint for clinical
        logger.info("\nAnnotating coding effects of mutations with varmint (clinical)")
        section_start = time.perf_counter()
        try:
            varmint(bam_files, fna_path, gff_path, dirs["tsv_clinical"], workers=cpu_count if args.fast else 4)
        except Exception as e:
            return f"Error running varmint on the alignment files: {e}"  
        logger.info(f"Varmint (clinical): {time.perf_counter() - section_start:.2f}s")

    # delete alignment files if not requested
    if not args.all:
        section_start = time.perf_counter()
        shutil.rmtree(dirs["alignment_dir"])
        shutil.rmtree(dirs["reference_dir"])
        logger.info(f"\nRemoved intermediate files")
        logger.info(f"Deleting intermediate files: {time.perf_counter() - section_start:.2f}s")

    # print run time
    cli_end_time = time.perf_counter()
    time_taken = round((cli_end_time - cli_start_time), 2) 
    logger.info(f"Run time of mutatea_cli: {timedelta(seconds=time_taken)}")

if __name__ == "__main__":
    mutatea()