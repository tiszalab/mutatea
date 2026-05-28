import gzip                                     # needed for unzipping reference files
import shutil                                   # needed for copying files
from Bio import SeqIO                           # needed for parsing fasta files
import glob                                     # needed for finding files
import os                                       # needed for file operations
import re                                       # needed for regular expressions (used for pulling out poolIDs)
import subprocess                               # needed for running shell commands
import json                                     # needed for parsing json files (custom dictionaries)
from multiprocessing import Pool                # needed for parallel processing
from variant_funcs import met_variant_alleles   # needed for variant labelling
import pysam                                    # needed for alingment quality filtering
import pandas as pd                             # needed for metadata processing

# process reference files
def process_reference_files(input_folder: str, reference_dir: str) -> tuple[str,str]:
    # make sure reference_dir exists
    os.makedirs(reference_dir, exist_ok=True)
    output_paths=[]

    # set empty strings to catch output
    fna_path = None
    gff_path = None

    # find all .fna and .gff files in the input folder
    unzipped_files = (
        glob.glob(os.path.join(input_folder, "*.fna"))
        + glob.glob(os.path.join(input_folder, "*.gff"))
    )

    # error for if there are multiple input files
    if len(unzipped_files) > 2:
        raise ValueError (f"only one gff(.gz) and one fna(.gz) can exist in the input folder")
    
    # copy any uncompressed reference files into the reference_dir
    for src_path in unzipped_files:
        filename = os.path.basename(src_path)
        out_path = os.path.join(reference_dir, filename)

        # make sure the reference file is not already in the reference_dir
        if os.path.exists(out_path):
            print(f"Existing file found in reference directory: {filename}")
            output_paths.append(out_path)
            continue
        if src_path != out_path:
            shutil.copy(src_path, out_path)
        print(f"Unzipped reference file was copied to: {out_path}")
        output_paths.append(out_path)

    # find all zipped files
    zipped_files = (
    glob.glob(os.path.join(input_folder, "*.fna.gz"))
    + glob.glob(os.path.join(input_folder, "*.gff.gz"))
    )

    # error for if there are multiple input files
    if len(zipped_files) > 2:
        raise ValueError (f"only one gff(.gz) and one fna(.gz) can exist in the input folder")  

    # process zipped reference files
    for zipped_file in zipped_files:
        filename = os.path.basename(zipped_file)[:-3]
        out_path = os.path.join(reference_dir, filename)

        # check if the unzipped version already exists in the reference dir
        if os.path.exists(out_path):
            output_paths.append(out_path)
            continue
        with gzip.open(zipped_file, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped reference file to: {out_path}")
        output_paths.append(out_path)

    # save paths to unzipped reference files
    fna_path = glob.glob(os.path.join(reference_dir, "*.fna"))[0]
    gff_path = glob.glob(os.path.join(reference_dir, "*.gff"))[0]

    # detailed error messages
    if not fna_path:
        raise ValueError(f"No fna(.gz) file found in the reference directory {reference_dir}")
    if not gff_path:
        raise ValueError(f"No gff(.gz) file found in the reference directory {reference_dir}")

    return fna_path, gff_path

# load in and merge metadata files
def process_metadata(metadata_folder:str, grouping:str = "month", logger=None) -> pd.DataFrame:
    metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))
    if not metadata_files:
        return pd.DataFrame()

    # filter to require specific columns
    required_columns = ["City", "Sample_ID", "Date"]
    md_list = []

    for file in metadata_files:
        df = pd.read_excel(file)
        # make sure all required columns are there
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"File {os.path.basename(file)} is missing required column: {col}")
        # save the files that passed filter to a list
        md_list.append(df)
    
    # merge into metadata
    metadata=pd.concat(md_list, ignore_index=True)

    # add column for time unit to metadata
    metadata["Date"] = pd.to_datetime(metadata["Date"], errors="coerce")

    # raise warning for rows with unparseable date formats
    bad_dates = metadata[metadata["Date"].isna()]
    if not bad_dates.empty:
        if logger: logger.warning(f"{len(bad_dates)} rows dropped due to unparseable Date")

    # crm: maybe need to add a filter, if there is no grouping column it shouldn't assume month_year (filter everywhere)
    if grouping == "day":
        metadata["Day_Month_Year"] = metadata["Date"].dt.strftime("%d_%m_%Y")
    elif grouping == "week":
        metadata["Week_Year"] = metadata["Date"].dt.strftime("%U_%Y")
    elif grouping == "year":
        metadata["Year"] = metadata["Date"].dt.strftime("%Y")
    else:
        metadata["Month_Year"] = metadata["Date"].dt.strftime("%m_%Y")
    
    # crm: not relevant to general user
    # add a sitecode column to metadata if not already present (older metadata files don't have this column)
    if "SiteCode" not in metadata.columns:
        metadata["SiteCode"] = pd.NA
    return metadata

# if include region: add region column to merged metadata
def add_region(metadata: pd.DataFrame, region_map_file: str = None) -> pd.DataFrame:
    # added filter for if there is not yet a city_region dictionary
    default_city_region = {
        "Houston, TX": "6_5S",
        "El Paso, TX": "9_10",
        "Lubbock, TX": "1",
        "Brownsville, TX": "11",
        "Wichita Falls, TX": "2_3",
        "Baytown, TX": "6_5S",
        "Humble, TX": "6_5S",
        "Missouri City, TX": "6_5S",
        "Austin, TX": "7",
        "Laredo, TX": "11",
        "Waco, TX": "7",
        "Fort Worth, TX": "2_3",
        "Palestine, TX": "4_5N",
        "Athens, TX": "4_5N",
        "Dallas, TX": "2_3",
        "DFW Airport, TX": "DFW_Airport",
        "Katy, TX": "6_5S",
        "San Antonio, TX": "8"
    }

    # use custom mapping if user inputted
    if region_map_file and os.path.exists(region_map_file):
        with open(region_map_file, 'r') as f:
            city_region = json.load(f)

        # replace spaces with underscores in region names (some counties are multiple words so they have spaces)
        city_region = {city: region.replace(' ', '_') for city, region in city_region.items()}

    else:
        city_region = default_city_region

    # use dictionary to add a "Region" column to the metadata
    metadata["Region"] = metadata["City"].map(city_region)
    # add an error for any unexpected cities, that way the user will know if they need to update the city_region dictionary
    unknown_cities = metadata.loc[metadata["Region"].isna(), "City"].unique()
    if len(unknown_cities) > 0:
        raise ValueError(f"Unknown cities found: {unknown_cities}. Please update the dictionary to include these cities.")
    else:
        # crm: need to adjust this, they may not be assigning cities to regions
        print("All cities in the metadata were successfully assigned to regions!\n")
    return metadata

# if include clinical: load in clinical metadata and fasta
def load_clinical_files(clinical_file_path: str, grouping:str = "month", logger=None) -> tuple[pd.DataFrame, str]:
    # find the csv in the clinical_metadata_path
    csv_file = glob.glob(os.path.join(clinical_file_path, "*.csv"))
    
    # raise error if no CSV files were found
    if not csv_file:
        raise FileNotFoundError(f"No CSV files found in {clinical_file_path}")
    # raise error if multiple CSV files were found
    if len(csv_file)>1:
        raise ValueError(f"Multiple CSV files found in {clinical_file_path}") 
    
    # read in csv
    # crm: added in low_memory to avoid dtypewarning
    clinical_metadata = pd.read_csv(csv_file[0], low_memory=False) 

    # set required columns for metadata csv
    required_columns = ["Accession", "Collection_Date"]

    # check for required columns
    for col in required_columns:
        if col not in clinical_metadata.columns:
            raise ValueError(f"Clinical csv is missing required column: {col}")

    # make sure the collection date column is a datetime object
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    
    # crm: adding fallback filter for month-based clinical data
    if grouping == "month":
        # keep clinical data with YYYY-MM values (no day) by converting them to the 1st of the month
        yr_mo_mask = clinical_metadata["Collection_Date"].isna()
        fallback = pd.to_datetime(
            clinical_metadata.loc[yr_mo_mask, "Collection_Date"].astype(str).str.strip() + "-01",
            format="%Y-%m-%d", errors="coerce"
        )
        clinical_metadata.loc[yr_mo_mask, "Collection_Date"] = fallback
        recovered = yr_mo_mask.sum() - clinical_metadata["Collection_Date"].isna().sum()
        
        # let user know we converted
        print(f"Recovered {recovered} rows with YYYY-MM dates (converted to 1st of month)")

    # raise warning for rows with unparseable date formats
    bad_dates = clinical_metadata[clinical_metadata["Collection_Date"].isna()]
    if not bad_dates.empty:
        if logger: logger.warning(f"{len(bad_dates)} rows dropped due to unparseable Collection_Date")
        if logger: logger.debug(f"Dropped Collection_Date value counts:\n{bad_dates['Collection_Date'].value_counts(dropna=False).to_string()}")

    # add unit of time column to the clinical metadata
    if grouping == "day":
        clinical_metadata["Day_Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%d_%m_%Y")
    elif grouping == "week":
        clinical_metadata["Week_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%U_%Y")
    elif grouping == "year":
        clinical_metadata["Year"] = clinical_metadata["Collection_Date"].dt.strftime("%Y")
    else:
        clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m_%Y")

    # find the fasta in the clinical_metadata_path
    fasta_file = glob.glob(os.path.join(clinical_file_path, "*.fasta"))
    # raise error if no fasta files were found
    if not fasta_file:
        raise FileNotFoundError(f"No fasta files found in {clinical_file_path}")
    # raise error if multiple fasta files were found
    if len(fasta_file)>1:
        raise ValueError(f"Multiple fasta files found in {clinical_file_path}")
    return clinical_metadata, fasta_file[0]

# create lists of accessions grouped by chosen unit of time
def create_grouped_accession_lists(clinical_metadata: pd.DataFrame, output_dir: str) -> None:
    # map grouping types to their corresponding column names
    grouping_columns = {
        "day": "Day_Month_Year",
        "week": "Week_Year", 
        "year": "Year",
        "month": "Month_Year" 
    }
    
    # set default of "Month_Year"
    # crm: do I need to set default here again?
    group_column = "Month_Year"
    
    # get column name by checking which grouping column exists in metadata
    for grouping_type, column_name in grouping_columns.items():
        if column_name in clinical_metadata.columns:
            group_column = column_name
            break
    
    # single loop for all grouping types
    for time, group in clinical_metadata.groupby(group_column):
        out_path = os.path.join(output_dir, f"{time}_list.txt")
        group["Accession"].to_csv(out_path, index=False, header=False)

# if include clinical: split clinical FASTA file by unit of time
def split_clinical_fasta_by_time(clinical_fasta_path: str, lists_dir: str, output_dir: str, logger=None) -> None:
    # load clinical fasta as dictionary
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))
    # build a stripped-key lookup once (e.g. PZ406333 -> PZ406333.1) for version-suffix mismatches
    stripped_to_full = {k.split(".")[0]: k for k in records_by_id}
    # loop through lists and split the clinical fasta by selected unit of time
    for list_file in glob.glob(os.path.join(lists_dir, "*.txt")):
        with open(list_file) as f:
            accessions = f.read().splitlines()

        # get all accessions for each unit of time
        # try exact match first, then fall back to stripped version suffix lookup
        time_accessions = []
        for a in accessions:
            if a in records_by_id:
                time_accessions.append(records_by_id[a])
            elif a in stripped_to_full:
                time_accessions.append(records_by_id[stripped_to_full[a]])

        # get the unit of time from the file name
        time = os.path.basename(list_file).replace("_list.txt", "")

        # skip if no accessions matched the fasta records
        if not time_accessions:
            if logger: logger.debug(f"Skipping {time}: no matching sequences found in clinical FASTA")
            continue

        # export clinical fasta by unit of time
        clinical_fasta_time = os.path.join(output_dir, f"{time}.fasta")
        SeqIO.write(time_accessions, clinical_fasta_time, "fasta")

# find wastewater reads from pools for the pathogen of interest
def find_wastewater_reads(pools_base_dir: str, pathogen: str, single_reads: bool = True) -> dict:
    # create empty dictionary to store reads by pool
    reads_by_pool = {}

    # for single reads
    if single_reads:
        # find all fasta and fastq files matching the pathogen
        fasta_files = sorted(glob.glob(os.path.join(pools_base_dir, f"*.{pathogen}.fasta")))
        fastq_files = sorted(glob.glob(os.path.join(pools_base_dir, f"*.{pathogen}.fastq")))
        fastq_gz_files = sorted(glob.glob(os.path.join(pools_base_dir, f"*.{pathogen}.fastq.gz")))
        all_files = fasta_files + fastq_files + fastq_gz_files
        
        if not all_files:
            print(f"No FASTA, FASTQ, or FASTQ.GZ files found for {pathogen} in {pools_base_dir}")
            return reads_by_pool
        
        # group files by pool_id (extracted from filename)
        for all_file in all_files:
            filename = os.path.basename(all_file)

            # crm: this could be an issue if they have these single reads in folders by poolID (not in the name of the fastq)
            # crm: may need to adjust for files that are not in folders with the p####

            # extract pool_id from filename
            parts = filename.split(".")
            pool_id = None
            for part in parts:
                if re.match(r'^p\d{4}$', part):
                    pool_id = part
                    break
            # crm: only processes files that have a valid pool_id
            if pool_id:            
                if pool_id not in reads_by_pool:
                    reads_by_pool[pool_id] = []
                reads_by_pool[pool_id].append(all_file)
    # for paired reads            
    else:
        known_r1_pattern = None
        known_r2_swap = None

        for pool_dir in sorted(glob.glob(os.path.join(pools_base_dir, "*"))):
            pool_id = os.path.basename(pool_dir)
            if not os.path.isdir(pool_dir) or not re.match(r'^p\d{4}$', pool_id):
                continue

            r1_files = []

            # if we already know the pattern, use it directly
            if known_r1_pattern:
                r1_files = sorted(glob.glob(os.path.join(pool_dir, "**", known_r1_pattern), recursive=True))
            else:
                # try each pattern until one hits
                r1_patterns = [
                    (f"*{pathogen}.R1.fastq",    (".R1.", ".R2.")),
                    (f"*{pathogen}.R1.fastq.gz", (".R1.", ".R2.")),
                    (f"*{pathogen}.R1.fasta",    (".R1.", ".R2.")),
                    (f"*{pathogen}_1.fastq",     ("_1.fastq", "_2.fastq")),
                    (f"*{pathogen}_1.fastq.gz",  ("_1.fastq.gz", "_2.fastq.gz")),
                ]
                for pattern, swap in r1_patterns:
                    hits = glob.glob(os.path.join(pool_dir, "**", pattern), recursive=True)
                    if hits:
                        r1_files = sorted(hits)
                        known_r1_pattern = pattern
                        known_r2_swap = swap
                        print(f"Detected read pattern: {pattern}")
                        break

            read_pairs = []
            if r1_files and known_r2_swap:
                for r1_file in r1_files:
                    r2_file = r1_file.replace(known_r2_swap[0], known_r2_swap[1])
                    if os.path.exists(r2_file):
                        read_pairs.append((r1_file, r2_file))
                    else:
                        print(f"No R2 file found for {r1_file}")

                # add them to the dictionary
                if read_pairs:
                    reads_by_pool[pool_id] = read_pairs

        # let user know if no reads were found for that pathogen in that pool
        if not reads_by_pool:
            print(f"No R1 files were found for {pathogen} in {pools_base_dir}")

    return reads_by_pool

# helper function for later alignment of wastewater reads
def _align_wastewater_reads(pool_id: str, read_files: list, fna_path: str, pools: str, pathogen: str, threads: int, minimap_preset: str = "sr", min_mapq: int = 0) -> list:
    # create list to capture output BAM file paths
    bam_files = []

    # create output directory for each pool
    pool_output_dir = os.path.join(pools, pool_id)
    os.makedirs(pool_output_dir, exist_ok=True)

    # align and sort wastewater reads
    for read_file in read_files:
        # paired reads
        if isinstance(read_file, tuple):
            r1_file, r2_file = read_file

            # crm: only works if the first item is the sampleid, need to confirm naming
            # get sample_ID
            filename = os.path.basename(r1_file)
            sample_name = filename.split(".")[0]

            minimap_cmd = ["minimap2", "-t", str(threads), "-ax", minimap_preset, fna_path, r1_file, r2_file]
        # single reads
        else:
            # extract sample_name from filename
            filename = os.path.basename(read_file)
            parts = filename.split(".")

            # remove file extension
            parts = parts[:-1]

            # remove pool ID if it's there
            if pool_id in parts:
                parts.remove(pool_id)

            # remove pathogen if it's there
            if pathogen:
                parts = [p for p in parts if p.lower() != pathogen.lower()]

            # get sample name
            sample_name = ".".join(parts) if parts else "unknown"

            minimap_cmd = ["minimap2", "-t", str(threads), "-ax", minimap_preset, fna_path, read_file]

        # create output BAM filename
        if min_mapq > 0:
            output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.mapq.sort.bam")
        else:
            output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")

        if os.path.exists(output_bam):
            with pysam.AlignmentFile(output_bam, "rb") as bam_cached:
                kept = bam_cached.count(until_eof=True)
            if kept > 0:
                bam_files.append(output_bam)
            continue

        try:
            # pipe minimap2 stdout into pysam, filter by mapq, write sorted BAM
            with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as mini2_proc:
                with pysam.AlignmentFile(mini2_proc.stdout, "r") as sam_stream:
                    unsorted_bam = output_bam.replace(".sort.bam", ".bam")
                    kept = 0
                    with pysam.AlignmentFile(unsorted_bam, "wb", header=sam_stream.header) as bam_out:
                        for read in sam_stream:
                            if read.mapping_quality >= min_mapq:
                                kept += 1
                                bam_out.write(read)

            if kept > 0:
                pysam.sort("-@", str(threads), "-o", output_bam, unsorted_bam)
                pysam.index(output_bam)
                bam_files.append(output_bam)

            os.remove(unsorted_bam)

        except Exception as e:
            print(f"Error processing {sample_name}: {e}")
            continue

    return bam_files

# align wastewater reads to reference files
def align_wastewater_reads(reads_by_pool: dict, fna_path: str, pools: str, pathogen: str, minimap_preset: str = "sr", threads: int = 8, workers: int = 4, min_mapq: int = 0) -> list:
    # create list to capture output
    bam_files = []
    
    if not reads_by_pool:
        return bam_files

    # prepare tasks
    tasks = []
    for pool_id, read_files in reads_by_pool.items():
        tasks.append((pool_id, read_files, fna_path, pools, pathogen, threads, minimap_preset, min_mapq))

    print(f"Aligning wastewater reads from {len(tasks)} pools using {workers} parallel workers")

    # run multiprocess
    with Pool(processes=workers) as pool:
        results = pool.starmap(_align_wastewater_reads, tasks)
    
    # combine BAM files by pool
    for pool_bam_files in results:
        bam_files.extend(pool_bam_files)
    
    return bam_files

# use wastewater metadata to group bam files by unit of time (option for if including region)
def create_wastewater_bam_groups(bam_files: list, metadata: pd.DataFrame, time_output_dir: str, region_output_dir: str = None, include_region: bool = True) -> str:
    # create empty dictionary to store file path lists
    bam_path_lists = {}
    
    # get bam directory from the inputted bam_files
    bam_dir = os.path.dirname(os.path.dirname(bam_files[0]))
    
    # build a sample_id -> bam_path lookup dictionary (once, instead of scanning per sample)
    sample_to_bam = {}
    for bam_file in bam_files:
        basename = os.path.basename(bam_file)
        # extract sample_id from basename (everything before the pool_id)
        sample_id = basename.split(".")[0]
        sample_to_bam[sample_id] = bam_file

    # map grouping types to their corresponding column names
    grouping_columns = {
        "day": "Day_Month_Year",
        "week": "Week_Year", 
        "year": "Year",
        "month": "Month_Year" 
    }
    
    # get column name (default is Month_Year)
    group_column = "Month_Year"  # default
    
    for grouping_type, column_name in grouping_columns.items():
        if column_name in metadata.columns:
            group_column = column_name
            break

    for time, group in metadata.groupby(group_column):
        # create empty list for the bam paths
        bam_paths=[]

        # loop through each row in the group and get the sample_id
        for sample_id in group["Sample_ID"]:
            # look up the bam path from the dictionary
            bam_path = sample_to_bam.get(sample_id)
            
            # make sure the file exists before adding
            if bam_path and os.path.exists(bam_path):
                bam_paths.append(bam_path)

        # add the list to the dictionary
        bam_path_lists[time] = bam_paths

    # loop through the bam path lists and create the text files
    for time, bam_path_list in bam_path_lists.items():
        # only create the txt file if the list contains bam paths
        if bam_path_list:
            out_path = os.path.join(time_output_dir, f"{time}_list.txt")
            with open(out_path, "w") as f:
                f.write("\n".join(bam_path_list))

    # do the same if region is also included
    if include_region and region_output_dir:
        # create empty dictionary to store bam path lists by region
        bam_path_lists_region = {}
    
        # time_region: loop through the metadata and create bam path lists
        for (time, region), group in metadata.groupby([group_column, "Region"]):
            # create empty list for the bam paths
            bam_paths_region = []

            # loop through each row in the group and get the sample_id
            for sample_id in group["Sample_ID"]:
                # look up the bam path from the dictionary
                bam_path = sample_to_bam.get(sample_id)
                
                # make sure the file exists before adding
                if bam_path and os.path.exists(bam_path):
                    bam_paths_region.append(bam_path)

            # create a combination time_region
            time_region = f"{time}_{region}"

            # add the list to the dictionary
            bam_path_lists_region[time_region] = bam_paths_region
     
        # loop through the bam path lists and create the list files
        for time_region, bam_path_list in bam_path_lists_region.items():
            # only create the txt file if the list contains bam paths
            if bam_path_list:
                out_path = os.path.join(region_output_dir, f"{time_region}_list.txt")
                with open(out_path, "w") as f:
                    f.write("\n".join(bam_path_list))
    if include_region:
        return time_output_dir, region_output_dir
    else:
        return time_output_dir

# merge bam files by time and time_region
def merge_wastewater_bams(list_dir: str, output_dir: str, threads: int = 8, min_mapq: int = 0) -> list:
    # create empty list to catch output
    merged_bams = []

    for list_file in glob.glob(os.path.join(list_dir, "*.txt")):
        # get the base name by removing the extension (can be for either time or time_region)
        list_name = os.path.basename(list_file).replace("_list.txt", "")

        # read BAM file paths from the list
        with open(list_file, "r") as f:
            bam_paths = f.read().splitlines()

        # create filename for outputted merged bam
        if min_mapq > 0:
            output_bam = os.path.join(output_dir, f"{list_name}.mapq.sort.bam")
        else:
            output_bam = os.path.join(output_dir, f"{list_name}.sort.bam")

        # samtools merge | samtools sort using list file to avoid argument list too long
        cmd = f"samtools merge -@ {threads} -f -b {list_file} - | samtools sort -@ {threads} -o {output_bam}"

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
                
            # index the sorted bam file
            subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)

            # only keep merged BAM if it contains aligned reads
            with pysam.AlignmentFile(output_bam, "rb") as bam:
                aligned = bam.count(until_eof=False)
            if aligned > 0:
                merged_bams.append(output_bam)
            else:
                os.remove(output_bam)
                bai = output_bam + ".bai"
                if os.path.exists(bai):
                    os.remove(bai)
            
        except subprocess.CalledProcessError as e:
            print(f"Error processing {list_name}: {e}")
            continue
    return merged_bams

## if include clinical: align clinical reads to reference
# helper function for later alignment of clinical reads
def _align_clinical_reads(fasta_file, fna_path, output_dir, threads, grouping, minimap_preset: str = "asm10", min_mapq: int = 0):

    # get the base name by removing the extension (can be for either time or time_region)
    time = os.path.basename(fasta_file).replace(".fasta", "")

    # catch output bam
    if min_mapq > 0:
        output_bam = os.path.join(output_dir, f"{time}.mapq.sort.bam")
    else:
        output_bam = os.path.join(output_dir, f"{time}.sort.bam")

    if os.path.exists(output_bam):
        with pysam.AlignmentFile(output_bam, "rb") as bam_cached:
            kept = bam_cached.count(until_eof=True)
        return output_bam if kept > 0 else None

    try:
        minimap_cmd = ["minimap2", "-ax", minimap_preset, fna_path, fasta_file]

        # pipe minimap2 stdout into pysam, filter by mapq, write sorted BAM
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as mini2_proc:
            with pysam.AlignmentFile(mini2_proc.stdout, "r") as sam_stream:
                unsorted_bam = output_bam.replace(".sort.bam", ".bam")
                kept = 0
                with pysam.AlignmentFile(unsorted_bam, "wb", header=sam_stream.header) as bam_out:
                    for read in sam_stream:
                        if read.mapping_quality >= min_mapq:
                            kept += 1
                            bam_out.write(read)

        if kept > 0:
            pysam.sort("-@", str(threads), "-o", output_bam, unsorted_bam)
            pysam.index(output_bam)

        os.remove(unsorted_bam)

        return output_bam if kept > 0 else None

    except Exception as e:
        print(f"Error processing {time}: {e}")
        return None

# align clinical reads to reference files
def align_clinical_reads(clinical_fasta_time:str, fna_path:str, output_dir: str, minimap_preset: str = "asm10", threads: int = 8, workers: int = 4, grouping: str = "month", min_mapq: int = 0) -> list:
    # create empty list to catch outputted bam files
    bam_files = []

    # find fasta files in the clinical fasta folder
    fasta_files = sorted(glob.glob(os.path.join(clinical_fasta_time, "*.fasta")))

    # prepare tasks
    tasks = []
    for fasta_file in fasta_files:
        # for all clinical reads, append the arguments to the tasks
        tasks.append((fasta_file, fna_path, output_dir, threads, grouping, minimap_preset, min_mapq))

    # print line is now saying number of tasks run with number of workers
    print(f"Aligning {len(tasks)} clinical fasta files using {workers} parallel workers")

    # run multiprocess 
    with Pool(processes=workers) as pool:
        results = pool.starmap(_align_clinical_reads, tasks)

    # combine BAM files
    for result in results:
        if result is not None:
            bam_files.append(result)
    
    # show any errors
    for result in results:
        if result.startswith("Error"):
            print(result)

    return bam_files

# optional statistics to get depth and breadth of genome coverage
def run_stats(bam_files:list, output_dir:str, logger=None) -> list:
    stats_files = []

    for bam_file in bam_files:
        # get base name from BAM file
        merge_name = os.path.basename(bam_file).replace(".mapq.sort.bam", "").replace(".sort.bam", "")

        try: 
            # filter: use number of reads aligned from each bam to only save stats files with content
            cmd_filter = ["samtools", "view", "-c", "-F", "4", bam_file]
            result = subprocess.run(cmd_filter, check=True, capture_output=True)
            aligned_reads = int(result.stdout.strip()) 

            # crm test to confirm that the stats file has content before saving
            if aligned_reads > 0:
                # create coverage file
                output_cov = os.path.join(output_dir, f"{merge_name}coverage.out")
                cmd_cov = ["samtools", "coverage", bam_file, "-o", output_cov]
                subprocess.run(cmd_cov, check=True, capture_output=True)
                stats_files.append(output_cov)
                
                # create stats file
                output_stats = os.path.join(output_dir, f"{merge_name}stats.out")
                cmd_stats = f"samtools stats {bam_file} > {output_stats}"
                subprocess.run(cmd_stats, shell=True, check=True)
                stats_files.append(output_stats)
            else:
                if logger: logger.debug(f"Skipping {merge_name}: contained no aligned reads")
                continue
                
        # error
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools on {merge_name}: {e}")
            # crm test
            continue
    return stats_files

# helper function for later running of varmint on merged bam files
def _varmint(bam_file, fna_path, gff_path, output_dir):
    # get the base name by removing the extension (can be for either time or time+region)
    merge_name = os.path.basename(bam_file).replace(".mapq.sort.bam", "").replace(".sort.bam", "")
    # create filename for outputted tsv
    output_tsv = os.path.join(output_dir, f"{merge_name}.tsv")

    # varmint
    try:
        df = met_variant_alleles(
            bam_path = bam_file,
            fasta_path = fna_path,
            gff_path = gff_path,
            min_base_qual = 20,
            min_depth = 1,
            min_map_qual =0
        )

        # confirm that the tsv file has content before saving
        if len(df) == 0:
            return f"Warning: {merge_name} produced empty TSV (deleted)"
        else:
            # write to TSV
            df.write_csv(output_tsv, separator="\t")
            return f"Success: {merge_name}"
    except Exception as e:
        return f"Error processing {merge_name}: {e}"

# run varmint on merged bam files
def varmint(bam_files: list, fna_path: str, gff_path: str, output_dir: str, workers: int = 4) -> None:    
    # prepare tasks
    tasks = []
    for bam_file in bam_files:
        # for all reads, append the arguments to the tasks
        tasks.append((bam_file, fna_path, gff_path, output_dir))
    
    # print line is now saying number of tasks run with number of workers, not number of reads/total per pool
    print(f"Running varmint using {workers} parallel workers")
    
    # run multiprocess 
    with Pool(processes=workers) as pool:
        results = pool.starmap(_varmint, tasks)
    
    # show any errors
    for result in results:
        if result.startswith("Error"):
            print(result)

def print_mutatea_banner():
    return r"""
                       _             _                  
   _ __ ___    _   _  | |_    __ _  | |_    ___    __ _ 
  | '_ ` _ \  | | | | | __|  / _` | | __|  / _ \  / _` |
  | | | | | | | |_| | | |_  | (_| | | |_  |  __/ | (_| |
  |_| |_| |_|  \__,_|  \__|  \__,_|  \__|  \___|  \__,_|
                                                                                                                                                                                                                                                                                                                             
            <████████████████████████████████>
         <█████$$$@&MMMMMMMMMMMMMMMMMMW@$$$█████>
      <█████$@&&WWWWMMM#==>MMMMMMMMMWWWWW&&@$█████>
   <████$%%&&&&WWWWMMMM*===*###*===>MMWWWW&&&&%%$████>
 <████@%%%&&&&WWWWMMM###<~~=##>~~==#MMMWWWW&&&&%%%@████>
<███@@%%%%&&&WWWWMMMM####<~~**~~~*###MMMWWWW&&&%%%%@@███>
$███@@@%%%%&&&WWWMM>===~=>><~~~~<*####MMMMWWW&&&&%%%@@@███$
$███@@@%%%&&&&WWWMM*===~~~~+~+~+~~~~~==>MMWWW&&&&%%%@@@███$  <ttttttt>
$███@@@%%%%&&&WWWMMMM####*>~+=~<>>~~===<MMWWW&&&&%%%@@@███$ <ttttttttt>
$████@@%%%%&&&WWWWMMMM##*~~~<*~~<#####MMMWWWW&&&%%%%@@████$<ttttttttttt>
$tt████@%%%&&&&WWWWMMMM<=~~>##=~~=##MMMMWWWW&&&&%%%@████tttttttttttttttt>
$tttt████$%%&&&&WWWWMM#===*####===<MMMMWWWW&&&&%%$████tttttttttttttttttt>
$ttttttt█████$@&&WWWWWMMMMMMMMM*==#MMWWWWW&&@$█████ttttttttt/    ttttttt>
$█tttttttttt███████@&WMMMMMMMMMMMMMMMW&@███████ttttttttttt█      ttttttt>
$t██tttttttttttttttt███████████████████ttttttttttttttttt██$      ttttttt>
$ttt██ttttttttttttttttttttttttttttttttttttttttttttttt██ttt$      ttttttt> 
$tttttt██ttttttttttttttttttttttttttttttttttttttttt███ttttt$     ttttttt>
$ttttttttt█████tttttttttttttttttttttttttttttt████ttttttttt$    ttttttt>
$ttttttttttttttt████████████████████████████tttttttttttttt$   ttttttt>
$ttttttttttttttttttttttttttttttttttttttttttttttttttttttttt$ tttttttt>
$█tttttttttttttttttttttttttttttttttttttttttttttttttttttttt█ttttttt>
$t█tttttttttttttttttttttttttttttttttttttttttttttttttttttt█ttttttt>
$tt███tttttttttttttttttttttttttttttttttttttttttttttttt███tttttt>
$ttttt████tttttttttttttttttttttttttttttttttttttttt████tttttttt>
$ttttttttt██████ttttttttttttttttttttttttttttt█████tttttttttt>
$tttttttttttttt██████████████████████████████ttttttttttttt>
 <ttttttttttttttttttttttttt|CRM|tttttttttttttttttttttttt>
  <ttttttttttttttttttttttttttttttttttttttttttttttttttt>
    <ttttttttttttttttttttttttttttttttttttttttttttttt>
       <tttttttttttttttttttttttttttttttttttttttttt>
           <ttttttttttttttttttttttttttttttttttt>
               <tttttttt|Tisza Lab|tttttttt>
"""
