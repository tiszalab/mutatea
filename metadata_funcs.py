import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
import glob
import os
import re
import subprocess

# load in and merge metadata files
def load_merge_metadata(metadata_folder:str) -> pd.DataFrame:
    metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))
    if not metadata_files:
        return pd.DataFrame()
    md_list=[pd.read_excel(file) for file in metadata_files]
    metadata=pd.concat(md_list, ignore_index=True)
    return metadata

# add month_year column to metadata
def add_month_year(metadata:pd.DataFrame) -> pd.DataFrame:
    metadata["Date"] = pd.to_datetime(metadata["Date"], errors="coerce")
    metadata["Month_Year"] = metadata["Date"].dt.strftime("%m.%Y")
    return metadata

# optional: add public health region column to merged metadata
def add_region(metadata: pd.DataFrame, city_region: dict = None) -> pd.DataFrame:
    # added filter for if there is not yet a city_region dictionary
    if city_region is None:
        city_region = {
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
    # use dictionary to add a "Region" column to the metadata
    metadata["Region"] = metadata["City"].map(city_region)
    # add a warning for any unexpected cities, that way the user will know if they need to update the city_region dictionary
    if len(metadata.loc[metadata["Region"].isna(), "City"].unique()) > 0:
        print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())
    else:
        print("\nAll cities in the metadata were successfully assigned to public health regions!\n")
    return metadata

# add a sitecode column to metadata if not already present (older metadata files don't have this column)
def ensure_sitecode_column(metadata: pd.DataFrame) -> pd.DataFrame:
    if "SiteCode" not in metadata.columns:
        metadata["SiteCode"] = pd.NA
    return metadata

# reorganize metadata columns to have a default order
def reorganize_metadata_columns(metadata: pd.DataFrame, no_region: bool = False) -> pd.DataFrame:
    if no_region == False:
        columns = [
            "City", "Sample_ID", "Site", "Date", "Flow",
            "PoolID", "SiteCode", "Region", "Month_Year"
        ]
    else:
        columns = [
            "City", "Sample_ID", "Site", "Date", "Flow",
            "PoolID", "SiteCode", "Month_Year"
        ]
    return metadata[columns]

# get the time range from the metadata
def get_date_range(metadata: pd.DataFrame) -> tuple:
    earliest_date = metadata["Date"].min()
    latest_date = metadata["Date"].max()
    return earliest_date, latest_date

# export the metadata as a csv to output path
def export_metadata(metadata: pd.DataFrame, metadata_dir: str, sep: str = ",") -> None:
    metadata.to_csv(f"{metadata_dir}/metadata_wastewater_combined.csv", sep=sep, index=False)

# load in clinical metadata and add column for month_year
def load_clinical_metadata(clinical_file_path: str) -> pd.DataFrame:
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
    # make sure the collection date column is a datetime object
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    # add month_year column to the clinical metadata
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")
    return clinical_metadata

# load in clinical fasta
def load_clinical_fasta(clinical_file_path: str) -> list:
    # find the fasta in the clinical_metadata_path
    fasta_file = glob.glob(os.path.join(clinical_file_path, "*.fasta"))
    # raise error if no fasta files were found
    if not fasta_file:
        raise FileNotFoundError(f"No fasta files found in {clinical_file_path}")
    # raise error if multiple fasta files were found
    if len(fasta_file)>1:
        raise ValueError(f"Multiple fasta files found in {clinical_file_path}")
    return fasta_file[0]

# process reference files
def process_reference_file(input_folder: str, reference_dir: str) -> list:
    # make sure reference_dir exists
    os.makedirs(reference_dir, exist_ok=True)
    output_paths=[]

    # find all .fna and .gff files in the input folder
    unzipped_files = (
        glob.glob(os.path.join(input_folder, "*.fna"))
        + glob.glob(os.path.join(input_folder, "*.gff"))
    )
    
    # first try to copy uncompressed reference files into the reference_dir
    if unzipped_files:
        for src_path in unzipped_files:
            filename = os.path.basename(src_path)
            out_path = os.path.join(reference_dir, filename)

            # make sure the reference file is not already in the reference_dir
            # crm: print line is clunky
            if os.path.exists(out_path):
                print(f"Existing reference file found in reference directory: {filename}")
                continue
            if src_path != out_path:
                shutil.copy(src_path, out_path)
            print(f"Reference file was already unzipped, was copied to: {out_path}")
            output_paths.append(out_path)

    # if no uncompressed files found then try to unzip files
    else:
        # find all zipped files
        gz_files = (
        glob.glob(os.path.join(input_folder, "*.fna.gz"))
        + glob.glob(os.path.join(input_folder, "*.gff.gz"))
        )

        # process zipped reference files
        for gz_file in gz_files:
            filename = os.path.basename(gz_file)[:-3]
            out_path = os.path.join(reference_dir, filename)
            with gzip.open(gz_file, "rb") as f_in, open(out_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"Unzipped reference file to: {out_path}")
            output_paths.append(out_path)

    # crm: being too specific about spacing, I know my character flaws
    print("")
    return output_paths

# create lists of accessions grouped by month_year
def create_monthly_accession_lists(clinical_metadata: pd.DataFrame, output_dir: str) -> None:
    for month, group in clinical_metadata.groupby("Month_Year"):
        out_path = Path(output_dir) / f"{month}_list.txt"
        group["Accession"].to_csv(out_path, index=False, header=False)

# split clinical FASTA file by monthly lists
def split_clinical_fasta_by_month(clinical_fasta_path: str, lists_dir: str, output_dir: str) -> None:
    # load clinical fasta as dictionary
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))
    # loop through monthly lists and split the clinical fasta
    for list_file in Path(lists_dir).glob("*.txt"):
        with open(list_file) as f:
            accessions = f.read().splitlines()
        # get all accessions for each month
        month_accessions = [records_by_id[a] for a in accessions if a in records_by_id]
        # get the month_year from the file name
        month_year = list_file.name.split("_")[0]
        # export clinical fasta by month
        clinical_fasta_month = os.path.join(output_dir, f"{month_year}.fasta")
        SeqIO.write(month_accessions, clinical_fasta_month, "fasta")

# find wastewater reads from pools for the subtype of interest
def find_wastewater_reads(pools_base_dir: str, subtype: str) -> dict:
    # create empty dictionary to store reads by pool
    reads_by_pool = {}
    for pool_dir in glob.glob(os.path.join(pools_base_dir, "*")):
        pool_id = os.path.basename(pool_dir)

        # skip the folder if it is not a directory or doesn't match the naming of the pools
        if not os.path.isdir(pool_dir) or not re.match(r'^p\d{4}$', pool_id):
            continue
        r1_files = glob.glob(os.path.join(pool_dir, "**", f"*{subtype}.R1.fastq"), recursive=True)

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
            # prints how many read pairs were found for each pool
            if len(read_pairs)==1:
                print(f"Pool {pool_id} contained {len(read_pairs)} {subtype} read pair")
            else:
                print(f"Pool {pool_id} contained {len(read_pairs)} {subtype} read pairs")

    # let user know if no reads were found for that subtype in that pool
    if not reads_by_pool:
        print(f"No R1 files were found for {subtype} in {pool_id}")

    return reads_by_pool

# align wastewater reads to reference files
def align_wastewater_reads(reads_by_pool: dict, reference_dir: str, pools: str, threads: int = 4) -> dict:

    # extract reference fasta from reference_files
    reference_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))[0]
    reference_gff = glob.glob(os.path.join(reference_dir, "*.gff"))[0]

    # bam files by pool
    bam_files_by_pool = {}

    # loop through the reads_by_pool dictionary and align the reads to the reference
    for pool_id, read_pairs in reads_by_pool.items():

        print(f"Aligning reads from pool {pool_id}")
        # create output directory for each pool
        pool_output_dir = os.path.join(pools, pool_id)
        os.makedirs(pool_output_dir, exist_ok=True)

        for read_pair in read_pairs:
            r1_file, r2_file = read_pair

            # get sample name from the filename of R1
            sample_name= os.path.basename(r1_file).split(".")[0]

            # create output BAM filename 
            output_bam = os.path.join(pool_output_dir, f"{sample_name}.{pool_id}.sort.bam")

            # crm: need to add in a skip if the BAM is already created

            # minimap2 | samtools view | samtools sort
            cmd = f"minimap2 -ax sr {reference_fasta} {r1_file} {r2_file} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output_bam}"
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
                
                # index the sorted bam files
                subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)
            
            except subprocess.CalledProcessError as e:
                print(f"Error processing {sample_name}: {e}")
                continue

# crm: do I end as None or dict?
# use wastewater metadata to create lists of bam filepaths for each month (optionally: and region)
def create_wastewater_bam_lists(metadata: pd.DataFrame, bam_dir: str, month_output_dir: str, region_output_dir: str = None, include_region: bool = True) -> None:
    # create empty dictionary to store file path lists
    bam_path_lists = {}
    
    # month: loop through the metadata and create bam path lists
    for month, group in metadata.groupby("Month_Year"):
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
        bam_path_lists[month] = bam_paths

    # loop through the bam path lists and create the text files
    for month, bam_path_list in bam_path_lists.items():
        # only create the txt file if the list contains bam paths
        if bam_path_list:
            out_path = os.path.join(month_output_dir, f"{month}_list.txt")
            with open(out_path, "w") as f:
                f.write("\n".join(bam_path_list))

    # do the same if region is also included
    if include_region and region_output_dir:
        # create empty dictionary to store bam path lists by region
        bam_path_lists_region = {}
    
        # month_region: loop through the metadata and create bam path lists
        for (month_year, region), group in metadata.groupby(["Month_Year", "Region"]):
            # create empty list for the bam paths
            bam_paths_region = []

            # loop through each row in the group and get the sample_id and pool_id
            for i in range(len(group)):
                sample_id = group.iloc[i]["Sample_ID"]
                pool_id = group.iloc[i]["PoolID"]

                # find file path for each sample
                bam_path = os.path.join(bam_dir, pool_id, f"{sample_id}.{pool_id}.sort.bam")

                # make sure the file exists before adding
                if os.path.exists(bam_path):
                    bam_paths_region.append(bam_path)

            # create a combination month_region
            month_region = f"{month_year}_{region}"

            # add the list to the dictionary
            bam_path_lists_region[month_region] = bam_paths_region
     
        # loop through the bam path lists and create the list files
        for month_region, bam_path_list in bam_path_lists_region.items():
            # only create the txt file if the list contains bam paths
            if bam_path_list:
                out_path = os.path.join(region_output_dir, f"{month_region}_list.txt")
                with open(out_path, "w") as f:
                    f.write("\n".join(bam_path_list))

# merge bam files using month and month_region lists
def merge_wastewater_bams(list_dir: str, output_dir: str, threads: int = 4) -> None:
    for list_file in glob.glob(os.path.join(list_dir, "*.txt")):
        # get the base name by removing the extension (can be for either month or month_region)
        list_name = os.path.basename(list_file).replace("_list.txt", "")

        # read BAM file paths from the list
        with open(list_file, "r") as f:
            bam_paths = f.read().splitlines()

        # create filename for outputted merged bam
        output_bam = os.path.join(output_dir, f"{list_name}.sort.bam")

        # samtools merge | samtools sort
        bam_paths_str = " ".join(bam_paths)
        cmd = f"samtools merge -@ {threads} -f - {bam_paths_str} | samtools sort -@ {threads} -o {output_bam}"

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
                
            # index the sorted bam file
            subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)
            
        except subprocess.CalledProcessError as e:
            print(f"Error processing {list_name}: {e}")
            continue

# run varmint on merged bam files
# def varmint_wastewater()

# optional: if include clinical, then align fasta files to reference
## pipe minimap2 into samtools sort, then index
def align_clinical_reads(clinical_fasta_month: str, output_dir: str, reference_dir: str, threads: int = 4) -> dict:

    # extract reference fasta from reference_files
    reference_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))[0]
    reference_gff = glob.glob(os.path.join(reference_dir, "*.gff"))[0]

    # create empty dictionary to catch outputted bam files
    bam_files_month = {}

    # loop through the fasta files in the clinical fasta folder
    for fasta_file in glob.glob(os.path.join(clinical_fasta_month, "*.fasta")):

        # get the base name by removing the extension (can be for either month or month_region)
        month_year = os.path.basename(fasta_file).replace(".fasta", "")

        # crm ?
        output_bam = os.path.join(output_dir, f"{month_year}.sort.bam")

        # crm: this might actually kill the script if the output already exists, is this skip necessary?
        # crm: need to add in a skip if the BAM is already created
        if os.path.exists(output_bam):
            bam_files_month[month_year]
            continue

        # minimap2 | samtools view | samtools sort
        # crm: want to replace this print line with a progress bar
        #print(f"Aligning samples from {month_year} to the reference genome")

        cmd = f"minimap2 -ax sr {reference_fasta} {fasta_file} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output_bam}"
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            
            # index the sorted bam files
            subprocess.run(["samtools", "index", output_bam], check=True, capture_output=True)
            
            # add the bam file to the dictionary
            bam_files_month[month_year] = output_bam

        except subprocess.CalledProcessError as e:
            print(f"Error processing {sample_name}: {e}")
            continue

    return bam_files_month

# optional: if include clinical, then run varmint on merged clinical bam files
# def varmint()