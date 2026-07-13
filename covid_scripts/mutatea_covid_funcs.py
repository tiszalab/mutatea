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
import pysam                                    # needed for alignment quality filtering
import tempfile                                 # needed for temporary files during bz2 decompression
import pandas as pd                             # needed for metadata processing

# crm: different functions
# all alignment functions were removed 
# parse_covid_positive_samples turns the sample list into a dictionary, with 3 values associated with the key sample_id
# find_bam_files_from_covid_samples turns the file path of the demix.out into a bam
# next functions are normal, create_wastewater_bam_groups 

# crm: want to change
# want it to take bam input, instead of turning the file path of the demix.out into a bam

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
        
        # make sure date is in datetime format
        df["Date"] = pd.to_datetime(df["Date"])
        md_list.append(df)
    
    # combine all metadata files
    metadata = pd.concat(md_list, ignore_index=True)

    # add time grouping columns
    metadata['Year'] = metadata['Date'].dt.year.astype(str)
    metadata['Month'] = metadata['Date'].dt.month.astype(str)
    metadata['Day'] = metadata['Date'].dt.day.astype(str)
    metadata['Month_Year'] = metadata['Month'] + '_' + metadata['Year']
    metadata['Day_Month_Year'] = metadata['Day'] + '_' + metadata['Month_Year']
    metadata['Week'] = metadata['Date'].dt.isocalendar().week.astype(str)
    metadata['Week_Year'] = metadata['Week'] + '_' + metadata['Year']

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

    # if user provided a custom dictionary file
    if region_map_file:
        with open(region_map_file, 'r') as f:
            city_region = json.load(f)
        city_region = {city: region.replace(' ', '_') for city, region in city_region.items()}

    else:
        city_region = default_city_region

    # use dictionary to add a "Region" column to the metadata
    metadata["Region"] = metadata["City"].map(city_region)
    
    # warn about unmapped cities
    unmapped_cities = metadata[metadata["Region"].isna()]["City"].unique()
    if len(unmapped_cities) > 0:
        print(f"Warning: The following cities were not mapped to regions: {unmapped_cities}")

    return metadata

# if include clinical: load in clinical metadata and fasta
def load_clinical_files(clinical_file_path: str, grouping:str = "month", logger=None) -> tuple[pd.DataFrame, str]:
    # find the csv in the clinical_metadata_path
    csv_file = glob.glob(os.path.join(clinical_file_path, "*.csv"))
    
    if not csv_file:
        raise FileNotFoundError(f"No CSV file found in {clinical_file_path}")
    
    # load in clinical metadata
    clinical_metadata = pd.read_csv(csv_file[0])
    
    # make sure date is in datetime format
    clinical_metadata["Date"] = pd.to_datetime(clinical_metadata["Date"])
    
    # add time grouping columns
    clinical_metadata['Year'] = clinical_metadata['Date'].dt.year.astype(str)
    clinical_metadata['Month'] = clinical_metadata['Date'].dt.month.astype(str)
    clinical_metadata['Day'] = clinical_metadata['Date'].dt.day.astype(str)
    clinical_metadata['Month_Year'] = clinical_metadata['Month'] + '_' + clinical_metadata['Year']
    clinical_metadata['Day_Month_Year'] = clinical_metadata['Day'] + '_' + clinical_metadata['Month_Year']
    clinical_metadata['Week'] = clinical_metadata['Date'].dt.isocalendar().week.astype(str)
    clinical_metadata['Week_Year'] = clinical_metadata['Week'] + '_' + clinical_metadata['Year']

    # find the fasta in the clinical_metadata_path
    fasta_file = glob.glob(os.path.join(clinical_file_path, "*.fasta"))
    
    if not fasta_file:
        raise FileNotFoundError(f"No FASTA file found in {clinical_file_path}")

    return clinical_metadata, fasta_file[0]

# create lists of accessions grouped by chosen unit of time
def create_grouped_accession_lists(clinical_metadata: pd.DataFrame, output_dir: str) -> None:
    # map grouping types to their corresponding column names
    grouping_columns = {
        "day": "Day_Month_Year",
        "week": "Week_Year", 
        "month": "Month_Year" 
    }
    
    # set default grouping of "Month_Year"
    group_column = "Month_Year"
    
    # get column name by checking which grouping column exists in metadata
    for grouping_type, column_name in grouping_columns.items():
        if column_name in clinical_metadata.columns:
            group_column = column_name
            break
    
    # group by the chosen time unit
    grouped = clinical_metadata.groupby(group_column)
    
    # create a list file for each group
    for group_name, group in grouped:
        out_path = os.path.join(output_dir, f"{group_name}.txt")
        group["Accession"].to_csv(out_path, index=False, header=False)

# if include clinical: split clinical FASTA file by unit of time
def split_clinical_fasta_by_time(clinical_fasta_path: str, lists_dir: str, output_dir: str, logger=None) -> None:
    # load clinical fasta as dictionary
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))
    
    # get all list files
    list_files = glob.glob(os.path.join(lists_dir, "*.txt"))
    
    for list_file in list_files:
        # get the time group from the filename
        time_group = os.path.basename(list_file).replace(".txt", "")
        
        # read the list of accessions
        with open(list_file, 'r') as f:
            accessions = [line.strip() for line in f if line.strip()]
        
        # collect records for this time group
        time_accessions = []
        for accession in accessions:
            if accession in records_by_id:
                time_accessions.append(records_by_id[accession])
        
        # write the time-grouped fasta
        if time_accessions:
            clinical_fasta_time = os.path.join(output_dir, f"{time_group}.fasta")
            SeqIO.write(time_accessions, clinical_fasta_time, "fasta")

# get read paths from covid sample list
# returns a dictionary with sample_id as the key and the sample info as value
def parse_covid_positive_samples(covid_samples_file: str) -> dict:
    samples_info = {}
    
    if not os.path.exists(covid_samples_file):
        raise FileNotFoundError(f"Covid positive samples file not found: {covid_samples_file}")
    
    with open(covid_samples_file, 'r') as f:
        for line in f:
            line = line.strip()
            # skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            # parse line: sample_id pool_id coverage demix_file
            # crm: may need to adjust, the sample list the user inputs may not be the same as the format I'm using
            # crm: maybe glob to find the column that has all the slashes? would indicate it contains the file paths
            parts = line.split()
            if len(parts) >= 4:
                sample_id = parts[0]
                pool_id = parts[1]
                coverage = float(parts[2])
                demix_file = parts[3]
                
                samples_info[sample_id] = {
                    'pool_id': pool_id,
                    'coverage': coverage,
                    'demix_file': demix_file
                }
    
    return samples_info

# crm: there is a sort.bam and a sort.bam.bz2 for each file, should only accept one input
# find BAM files from covid positive samples
def find_bam_files_from_covid_samples(covid_samples_file: str, bam_base_dir: str, min_mapq: int = 0) -> list:
    # get covid positive samples info
    samples_info = parse_covid_positive_samples(covid_samples_file)
    
    # transform demix file paths to BAM file paths
    bam_files = []
    for sample_id, sample_info in samples_info.items():
        demix_file = sample_info['demix_file']
        # crm: try .sort.bam first, then .sort.bam.bz2
        bam_file = demix_file.replace('.demix.out', '.sort.bam')
        bam_file_bz2 = demix_file.replace('.demix.out', '.sort.bam.bz2')
        
        if os.path.exists(bam_file):
            bam_files.append(bam_file)
        elif os.path.exists(bam_file_bz2):
            bam_files.append(bam_file_bz2)
        else:
            print(f"Warning: BAM file not found for {sample_id}: {bam_file}")
    
    print(f"Found {len(bam_files)} BAM files from {len(samples_info)} COVID-positive samples")
    
    # filter BAM files by MAPQ if specified
    if min_mapq > 0:
        print(f"Filtering BAM files by MAPQ >= {min_mapq}")
        filtered_bams = []
        for bam_file in bam_files:
            filtered_bam = bam_file.replace(".bam", f".mapq{min_mapq}.bam")
            cmd_filter = f"samtools view -b -q {min_mapq} {bam_file} > {filtered_bam}"
            try:
                subprocess.run(cmd_filter, shell=True, check=True, executable="/bin/bash")
                filtered_bams.append(filtered_bam)
            except Exception as e:
                print(f"Warning: Could not filter {bam_file}: {e}")
                filtered_bams.append(bam_file)
        bam_files = filtered_bams
        print(f"Filtered to {len(bam_files)} BAM files")
    
    return bam_files

# use wastewater metadata to group bam files by unit of time (option for if including region)
def create_wastewater_bam_groups(bam_files: list, metadata: pd.DataFrame, time_output_dir: str, region_output_dir: str = None, include_region: bool = True, grouping: str = "month") -> str:
    # create empty dictionary to store file path lists
    bam_path_lists = {}

    # build a sample_id -> bam_path lookup dictionary keyed by sample_id from filename
    sample_to_bam = {}
    for bam_file in bam_files:
        basename = os.path.basename(bam_file)
        sample_id = basename.split(".")[0]
        sample_to_bam[sample_id] = bam_file

    # map grouping types to their corresponding column names
    grouping_columns = {
        "day": "Day_Month_Year",
        "week": "Week_Year",
        "year": "Year",
        "month": "Month_Year"
    }

    # use the specified grouping to select the correct column
    group_column = grouping_columns.get(grouping, "Month_Year")

    for time, group in metadata.groupby(group_column):
        bam_paths = []
        for sample_id in group["Sample_ID"]:
            bam_path = sample_to_bam.get(sample_id)
            if bam_path and os.path.exists(bam_path):
                bam_paths.append(bam_path)
        bam_path_lists[time] = bam_paths

    for time, bam_path_list in bam_path_lists.items():
        if bam_path_list:
            out_path = os.path.join(time_output_dir, f"{time}_list.txt")
            with open(out_path, "w") as f:
                f.write("\n".join(bam_path_list))

    if include_region and region_output_dir:
        bam_path_lists_region = {}
        for (time, region), group in metadata.groupby([group_column, "Region"]):
            bam_paths_region = []
            for sample_id in group["Sample_ID"]:
                bam_path = sample_to_bam.get(sample_id)
                if bam_path and os.path.exists(bam_path):
                    bam_paths_region.append(bam_path)
            time_region = f"{time}_{region}"
            bam_path_lists_region[time_region] = bam_paths_region

        for time_region, bam_path_list in bam_path_lists_region.items():
            if bam_path_list:
                out_path = os.path.join(region_output_dir, f"{time_region}_list.txt")
                with open(out_path, "w") as f:
                    f.write("\n".join(bam_path_list))

    if include_region:
        return time_output_dir, region_output_dir
    else:
        return time_output_dir

# merge bam files by time and time_region, unzip if they're sort.bam.bz2
# crm: need to confirm we're not just writing over the sort.bam when unzipping the sort.bam.bz2
def merge_wastewater_bams(list_dir: str, output_dir: str, threads: int = 8, min_mapq: int = 0) -> list:
    # create empty list to catch output
    merged_bams = []
    
    for list_file in glob.glob(os.path.join(list_dir, "*.txt")):
        group_name = os.path.basename(list_file).replace("_list.txt", "").replace(".txt", "")
        
        with open(list_file, 'r') as f:
            bam_files = [line.strip() for line in f if line.strip()]
        
        if not bam_files:
            continue
        
        temp_dir = tempfile.mkdtemp()
        try:
            # decompress any .bam.bz2 files to temp dir (samtools can't read bzip2 directly)
            ready_bams = []
            for bam_path in bam_files:
                if bam_path.endswith('.bam.bz2'):
                    temp_bam = os.path.join(temp_dir, os.path.basename(bam_path[:-4]))
                    subprocess.run(f"bzip2 -dkc {bam_path} > {temp_bam}", shell=True, check=True, executable="/bin/bash")
                    subprocess.run(["samtools", "index", temp_bam], check=True)
                    ready_bams.append(temp_bam)
                else:
                    ready_bams.append(bam_path)
            
            # write updated list with decompressed paths
            temp_list = os.path.join(temp_dir, "list.txt")
            with open(temp_list, 'w') as f:
                f.write("\n".join(ready_bams))
            
            # merge and sort and index
            output_bam = os.path.join(output_dir, f"{group_name}.sort.bam")
            cmd = f"samtools merge -@ {threads} -f -b {temp_list} - | samtools sort -@ {threads} -o {output_bam}"
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            subprocess.run(["samtools", "index", output_bam], check=True)
            
            # filter by MAPQ if requested
            if min_mapq > 0:
                filtered_bam = output_bam.replace(".sort.bam", f".mapq{min_mapq}.sort.bam")
                subprocess.run(f"samtools view -b -q {min_mapq} {output_bam} | samtools sort -@ {threads} -o {filtered_bam} -", shell=True, check=True, executable="/bin/bash")
                subprocess.run(["samtools", "index", filtered_bam], check=True)
                os.remove(output_bam)
                if os.path.exists(output_bam + ".bai"):
                    os.remove(output_bam + ".bai")
                merged_bams.append(filtered_bam)
            else:
                merged_bams.append(output_bam)
                
        except subprocess.CalledProcessError as e:
            print(f"Error merging BAM files for {group_name}: {e}")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    return merged_bams


# optional statistics to get depth and breadth of genome coverage
def run_stats(bam_files:list, output_dir:str, logger=None) -> list:
    stats_files = []

    for bam_file in bam_files:
        # get the base name by removing the extension
        base_name = os.path.basename(bam_file).replace(".sort.bam", "").replace(".mapq.sort.bam", "")
        # create filename for outputted tsv
        stats_file = os.path.join(output_dir, f"{base_name}_stats.tsv")
        
        # run samtools coverage
        cmd = f"samtools coverage {bam_file} > {stats_file}"
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            stats_files.append(stats_file)
        except subprocess.CalledProcessError as e:
            print(f"Error getting coverage statistics for {base_name}: {e}")

    return stats_files

# helper function for later running of varmint on merged bam files
def _varmint(bam_file, fna_path, gff_path, output_dir):
    # get the base name by removing the extension (can be for either time or time+region)
    merge_name = os.path.basename(bam_file).replace(".mapq.sort.bam", "").replace(".sort.bam", "")
    # create filename for outputted tsv
    output_tsv = os.path.join(output_dir, f"{merge_name}.tsv")

    try:
        df = met_variant_alleles(
            bam_path = bam_file,
            fasta_path = fna_path,
            gff_path = gff_path,
            min_map_qual = 0
        )

        if len(df) == 0:
            return f"Warning: {merge_name} produced empty TSV (deleted)"
        else:
            df.write_csv(output_tsv, separator="\t")
            return f"Success: {merge_name}"
    except Exception as e:
        return f"Error processing {merge_name}: {e}"

# run varmint on merged bam files
def varmint(bam_files: list, fna_path: str, gff_path: str, output_dir: str, workers: int = 4) -> None:    
    # prepare tasks
    tasks = []
    for bam_file in bam_files:
        tasks.append((bam_file, fna_path, gff_path, output_dir))
    
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
