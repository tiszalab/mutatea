import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
import glob
import os
import re
import subprocess
import json
from multiprocessing import Pool
from variant_funcs import met_variant_alleles

# process reference files
def process_reference_file(input_folder: str, reference_dir: str) -> tuple[str,str]:
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
def process_metadata(metadata_folder:str, grouping:str = "month") -> pd.DataFrame:
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

    # crm: maybe need to add a filter, if there is no grouping column it shouldn't assume month_year (filter everywhere)
    if grouping == "day":
        metadata["Day_Month_Year"] = metadata["Date"].dt.strftime("%d_%m_%Y")
    elif grouping == "week":
        metadata["Week_Year"] = metadata["Date"].dt.strftime("%U_%Y")
    elif grouping == "year":
        metadata["Year"] = metadata["Date"].dt.strftime("%Y")
    else:
        metadata["Month_Year"] = metadata["Date"].dt.strftime("%m_%Y")
    
    # add a sitecode column to metadata if not already present (older metadata files don't have this column)
    if "SiteCode" not in metadata.columns:
        metadata["SiteCode"] = pd.NA
    return metadata

# if include region: add public health region column to merged metadata
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
        "DFW Airport, TX": "2_3",
        "Katy, TX": "6_5S"
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
    # add a warning for any unexpected cities, that way the user will know if they need to update the city_region dictionary
    if len(metadata.loc[metadata["Region"].isna(), "City"].unique()) > 0:
        print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())
    else:
        print("All cities in the metadata were successfully assigned to public health regions!\n")
    return metadata

# if include clinical: load in clinical metadata and fasta
def load_clinical_files(clinical_file_path: str, grouping:str = "month") -> tuple[pd.DataFrame, str]:
    # find the csv in the clinical_metadata_path
    csv_file = glob.glob(os.path.join(clinical_file_path, "*.csv"))
    
    # raise error if no CSV files were found
    if not csv_file:
        raise FileNotFoundError(f"No CSV files found in {clinical_file_path}")
    # raise error if multiple CSV files were found
    if len(csv_file)>1:
        raise ValueError(f"Multiple CSV files found in {clinical_file_path}") 
    
    # read in csv
    clinical_metadata = pd.read_csv(csv_file[0]) 

    # set required columns for metadata csv
    required_columns = ["Accession", "Collection_Date"]

    # check for required columns
    for col in required_columns:
        if col not in clinical_metadata.columns:
            raise ValueError(f"Clinical csv is missing required column: {col}")

    # make sure the collection date column is a datetime object
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    
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
        out_path = Path(output_dir) / f"{time}_list.txt"
        group["Accession"].to_csv(out_path, index=False, header=False)

# if include clinical: split clinical FASTA file by unit of time
def split_clinical_fasta_by_time(clinical_fasta_path: str, lists_dir: str, output_dir: str) -> None:
    # load clinical fasta as dictionary
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))
    # loop through lists and split the clinical fasta by selected unit of time
    for list_file in Path(lists_dir).glob("*.txt"):
        with open(list_file) as f:
            accessions = f.read().splitlines()

        # get all accessions for each unit of time
        time_accessions = [records_by_id[a] for a in accessions if a in records_by_id]

        # get the unit of time from the file name
        time = list_file.name.split(".")[0]

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
        all_files = fasta_files + fastq_files
        
        if not all_files:
            print(f"No FASTA or FASTQ files found for {pathogen} in {pools_base_dir}")
            return reads_by_pool
        
        # group files by pool_id (extracted from filename)
        for all_file in all_files:
            filename = os.path.basename(all_file)

            # crm: this could be an issue if they have these single reads in folders by poolID (not in the name of the fastq)
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
        for pool_dir in sorted(glob.glob(os.path.join(pools_base_dir, "*"))):
            pool_id = os.path.basename(pool_dir)

            # crm: this may not be the way most people name their pools
            # skip the folder if it is not a directory or doesn't match the naming of the pools
            if not os.path.isdir(pool_dir) or not re.match(r'^p\d{4}$', pool_id):
                continue
            r1_files = glob.glob(os.path.join(pool_dir, "**", f"*{pathogen}.R1.fastq"), recursive=True)

            # create empty list for the paired reads
            read_pairs=[]
            # loop through all existing r1 reads to find their r2 read
            if r1_files:
                # crm: this is assuming that the respective r2 would be kept next to r1 (same pool), but maybe need to confirm that
                for r1_file in r1_files:
                    r2_file = r1_file.replace("R1.fastq", "R2.fastq")

                    # make sure r2 actually exists
                    if os.path.exists(r2_file):
                        read_pairs.append((r1_file, r2_file))
                    else:
                        print(f"No R2 file found for {r1_file}")

            # add them to the dictionary
            if read_pairs:
                reads_by_pool[pool_id] = read_pairs

        # let user know if no reads were found for that pathogen in that pool
        if not reads_by_pool:
            print(f"No R1 files were found for {pathogen} in {pool_id}")

    return reads_by_pool

# helper function for later alignment of wastewater reads
def _align_wastewater_reads(pool_id: str, read_files: list, fna_path: str, pools: str, pathogen: str, threads: int) -> list:
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

            # get sample name from the filename of R1
            filename = os.path.basename(r1_file)
            parts = filename.split(".")  

            # crm: only works if the first item is the sampleid, need to confirm naming
            # get sample_ID         
            sample_name = parts[0]

            # create output BAM filename 
            output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")

            # minimap2 | samtools view | samtools sort
            cmd = f"minimap2 -t {threads} -ax sr {fna_path} {r1_file} {r2_file} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output_bam}"  
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

            # create output BAM filename 
            output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")

            # minimap2 | samtools view | samtools sort
            cmd = f"minimap2 -t {threads} -ax sr {fna_path} {read_file} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output_bam}"
        
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
                                    
            # index the sorted bam files
            subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)
            
            # add to list of BAM files
            bam_files.append(output_bam)
                                    
        except subprocess.CalledProcessError as e:
            print(f"Error processing {sample_name}: {e}")
            continue
            
    return bam_files

# align wastewater reads to reference files
def align_wastewater_reads(reads_by_pool: dict, fna_path: str, pools: str, pathogen: str, threads: int = 8, workers: int = 4) -> list:
    # create list to capture output
    bam_files = []
    
    if not reads_by_pool:
        return bam_files

    # prepare tasks
    tasks = []
    for pool_id, read_files in reads_by_pool.items():
        tasks.append((pool_id, read_files, fna_path, pools, pathogen, threads))

    print(f"Aligning wastewater reads from {len(tasks)} pools using {workers} parallel workers")

    # run multiprocess
    with Pool(processes=workers) as pool:
        results = pool.starmap(_align_wastewater_reads, tasks)
    
    # combine BAM files by pool
    for pool_bam_files in results:
        bam_files.extend(pool_bam_files)
    
    return bam_files

# use wastewater metadata to group bam files by unit of time (option for if including region)
def create_wastewater_bam_groups(bam_files: list, metadata: pd.DataFrame, month_output_dir: str, region_output_dir: str = None, include_region: bool = True) -> str:
    # create empty dictionary to store file path lists
    bam_path_lists = {}
    
    # get bam directory from the inputted bam_files
    bam_dir = os.path.dirname(os.path.dirname(bam_files[0]))
    
    # time: loop through the metadata and create bam path lists
    # map grouping types to their corresponding column names
    grouping_columns = {
        "day": "Day_Month_Year",
        "week": "Week_Year", 
        "year": "Year",
        "month": "Month_Year" 
    }
    
    # get column name (default is Month_Year)
    group_column = "Month_Year"  # default
    max_unique = 0
    
    for grouping_type, column_name in grouping_columns.items():
        if column_name in metadata.columns:
            group_column = column_name
            break

    for time, group in metadata.groupby(group_column):
        # create empty list for the bam paths
        bam_paths=[]

        # loop through each row in the group and get the sample_id and pool_id
        for i in range(len(group)):
            sample_id = group.iloc[i]["Sample_ID"]
            pool_id = group.iloc[i]["PoolID"]

            # find file path for each sample
            bam_path = os.path.join(bam_dir, pool_id, f"{sample_id}.{pool_id}.sort.bam")

            # make sure the file exists before adding
            if os.path.exists(bam_path):
                bam_paths.append(bam_path)
        # add the list to the dictionary
        bam_path_lists[time] = bam_paths

    # loop through the bam path lists and create the text files
    for time, bam_path_list in bam_path_lists.items():
        # only create the txt file if the list contains bam paths
        if bam_path_list:
            out_path = os.path.join(month_output_dir, f"{time}_list.txt")
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

            # loop through each row in the group and get the sample_id and pool_id
            # crm: make sure you can explain iloc
            for i in range(len(group)):
                sample_id = group.iloc[i]["Sample_ID"]
                pool_id = group.iloc[i]["PoolID"]

                # find file path for each sample
                bam_path = os.path.join(bam_dir, pool_id, f"{sample_id}.{pool_id}.sort.bam")

                # make sure the file exists before adding
                if os.path.exists(bam_path):
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
        return month_output_dir, region_output_dir
    else:
        return month_output_dir

# merge bam files by time and time_region
def merge_wastewater_bams(list_dir: str, output_dir: str, threads: int = 8) -> list:
    # create empty list to catch output
    merged_bams = []

    for list_file in glob.glob(os.path.join(list_dir, "*.txt")):
        # get the base name by removing the extension (can be for either time or time_region)
        list_name = os.path.basename(list_file).replace("_list.txt", "")

        # read BAM file paths from the list
        with open(list_file, "r") as f:
            bam_paths = f.read().splitlines()

        # create filename for outputted merged bam
        output_bam = os.path.join(output_dir, f"{list_name}.sort.bam")

        # samtools merge | samtools sort
        # bam_paths_str = " ".join(bam_paths)
        cmd = f"samtools merge -@ {threads} -b - {list_file} | samtools sort -@ {threads} -o {output_bam}"

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
                
            # index the sorted bam file
            subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)

            # add to list of created bams
            merged_bams.append(output_bam)
            
        except subprocess.CalledProcessError as e:
            print(f"Error processing {list_name}: {e}")
            continue
    return merged_bams

## if include clinical: align clinical reads to reference
# helper function for later alignment of clinical reads
def _align_clinical_reads(fasta_file, fna_path, output_dir, threads, grouping):


    # get the base name by removing the extension (can be for either time or time_region)
    month_year = os.path.basename(fasta_file).replace(".fasta", "")

    # catch output bam
    output_bam = os.path.join(output_dir, f"{month_year}.sort.bam")

    # crm: do I need to add in a skip if the BAM is already created?
    if os.path.exists(output_bam):
        return output_bam
        
    try:
        # crm: changed this to asm 10, not sure how it will work
        # minimap2 | samtools view | samtools sort
        cmd = f"minimap2 -ax asm10 {fna_path} {fasta_file} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output_bam}"
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
            
        # index the sorted bam files
        subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)
            
        return output_bam
       
    # error
    except subprocess.CalledProcessError as e:
        print(f"Error processing {month_year}: {e}")
        return None

# align clinical reads to reference files
def align_clinical_reads(clinical_fasta_month:str, fna_path:str, output_dir: str, threads: int = 8, workers: int = 4, grouping: str = "month") -> list:
    # create empty list to catch outputted bam files
    bam_files = []

    # find fasta files in the clinical fasta folder
    fasta_files = sorted(glob.glob(os.path.join(clinical_fasta_month, "*.fasta")))

    # prepare tasks
    tasks = []
    for fasta_file in fasta_files:
        # for all clinical reads, append the arguments to the tasks
        tasks.append((fasta_file, fna_path, output_dir, threads, grouping))

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

# helper function for later running of varmint on merged bam files
def _varmint(bam_file, fna_path, gff_path, output_dir):
    # get the base name by removing the extension (can be for either month or month_region)
    merge_name = os.path.basename(bam_file).replace(".sort.bam", "")

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
    print(f"Running varmint on {len(tasks)} BAM files using {workers} parallel workers")
    
    # run multiprocess 
    with Pool(processes=workers) as pool:
        results = pool.starmap(_varmint, tasks)
    
    # show any errors
    for result in results:
        if result.startswith("Error"):
            print(result)