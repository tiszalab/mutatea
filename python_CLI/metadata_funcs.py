import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
from subprocess import Popen
import argparse
import glob
import os

# crm: need filter to only merge files that have specific columns
# load in and merge metadata files
def load_merge_metadata(metadata_folder:str) -> pd.DataFrame:
    metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))
    md_list=[pd.read_excel(file) for file in metadata_files]
    metadata=pd.concat(md_list, ignore_index=True)
    return metadata
        

# add month_year column to metadata
def add_month_year(metadata:pd.DataFrame) -> pd.DataFrame:
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
        print("All cities in the metadata were successfully assigned to public health regions!")
    return metadata

# add a sitecode column to metadata if not already present (older metadata files don't have this column)
def ensure_sitecode_column(metadata: pd.DataFrame) -> pd.DataFrame:
    if "SiteCode" not in metadata.columns:
        metadata["SiteCode"] = pd.NA
    return metadata

# crm: need to make sure this function works with default of region + month_year
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


# export the metadata as a csv to output path
def export_metadata(metadata: pd.DataFrame, metadata_dir: str, sep: str = ",") -> None:
    metadata.to_csv(f"{metadata_dir}/metadata_wastewater_combined.csv", sep=sep, index=False)

# get the time range from the metadata
def get_date_range(metadata: pd.DataFrame) -> tuple:
    earliest_date = metadata["Date"].min()
    latest_date = metadata["Date"].max()
    return earliest_date, latest_date

# load in clinical metadata and add column for month_year
def load_clinical_metadata(clinical_metadata_path: str) -> pd.DataFrame:
    clinical_metadata = pd.read_csv(clinical_metadata_path)
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")
    return clinical_metadata


###### crm: this is not functional, needs a lot of work
# process reference files
def process_reference_file(input_folder: str, reference_dir: str) -> list:
    # make sure reference_dir exists
    os.makedirs(reference_dir, exist_ok=True)
    output_paths=[]

    # unzip gz files into the reference_dir
    # find all .gz files in the input folder
    gz_files = glob.glob(os.path.join(input_folder, "*.gz"))
    for gz_file in gz_files:
        # Get the output filename by removing .gz extension
        filename = os.path.basename(gz_file)[:-3]
        output_path = os.path.join(reference_dir, filename)
        # unzip the file and copy it to the reference_dir
        with gzip.open(gz_file, "rb") as f_in, open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"\nUnzipped to: {output_path}\n")
        output_paths.append(output_path)
    # if the reference files are already unzipped, copy them to the reference_dir
    if ".gff" not in output_path:
        for pattern in ("*.fasta", "*.fna", "*.gff"):
            for src_path in glob.glob(os.path.join(input_folder, pattern)):
                filename = os.path.basename(src_path)
                output_path = os.path.join(reference_dir, filename)
                # make sure the reference file is not already in the reference_dir
                if os.path.exists(output_path):
                    print(f"Reference file {filename} already exists in {reference_dir}.\n")
                    continue
                if src_path != output_path:
                    shutil.copy(src_path, output_path)
                print(f"Copied to: {output_path}\n")
                output_paths.append(output_path)
    return output_paths

# split clinical fasta by month
def split_clinical_fasta(clinical_fasta:str, output_dir:str):
    pass


###################### AI versions, need to be vetted ######################

def copy_reference_file(input_path: str, output_dir: str) -> str:
    """Copy a reference file to the output directory if not already there."""
    filename = os.path.basename(input_path)
    output_path = os.path.join(output_dir, filename)
    if input_path != output_path:
        shutil.copy(input_path, output_path)
    return output_path


def process_reference_file(input_path: str, output_dir: str) -> str:
    """Process a reference file - unzip if needed, or copy to output dir."""
    if input_path.endswith(".gz"):
        return unzip_file(input_path, output_dir)
    else:
        return copy_reference_file(input_path, output_dir)


def create_monthly_accession_lists(clinical_metadata: pd.DataFrame, output_dir: str) -> None:
    """Create text files with accessions grouped by Month_Year."""
    for month, group in clinical_metadata.groupby("Month_Year"):
        out_path = Path(output_dir) / f"{month}_list.txt"
        group["Accession"].to_csv(out_path, index=False, header=False)


def split_clinical_fasta_by_month(clinical_fasta_path: str, lists_dir: str, output_dir: str) -> None:
    """Split a clinical FASTA file by monthly accession lists."""
    # Load clinical fasta as dictionary
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))
    
    # Loop through monthly lists and split the clinical fasta
    for list_file in Path(lists_dir).glob("*.txt"):
        with open(list_file) as f:
            accessions = f.read().splitlines()
        
        # Get all accessions for each month
        month_accessions = [records_by_id[a] for a in accessions if a in records_by_id]
        
        # Get the month_year from the file name
        month_year = list_file.name.split("_")[0]
        
        # Export clinical fasta by month
        clinical_monthly_fasta = os.path.join(output_dir, f"{month_year}.fasta")
        SeqIO.write(month_accessions, clinical_monthly_fasta, "fasta")

# create output directories
def create_output_directories(output_dir: str, include_clinical: bool = False) -> dict:
    dirs = {}
    
    # main output directory
    os.makedirs(output_dir, exist_ok=True)
    dirs["output"] = output_dir
    
    # cleaned and merged metadata directory
    metadata_dir = os.path.join(output_dir, "metadata_files")
    os.makedirs(metadata_dir, exist_ok=True)
    dirs["metadata"] = metadata_dir
    
    # output file directory
    reference_dir = os.path.join(output_dir, "reference_files")
    os.makedirs(reference_dir, exist_ok=True)
    dirs["reference"] = reference_dir
    
    # clinical directories
    if include_clinical:
        # parent folder for clinical output
        clinical_output = os.path.join(output_dir, "clinical_output")
        os.makedirs(clinical_output, exist_ok=True)
        dirs["clinical_output"] = clinical_output
        
        # folder for the lists of accessions by month
        clinical_lists = os.path.join(clinical_output, "monthly_lists")
        os.makedirs(clinical_lists, exist_ok=True)
        dirs["clinical_lists"] = clinical_lists
        
        # folder for the clinical fastas split by month
        clinical_fasta = os.path.join(clinical_output, "monthly_fasta")
        os.makedirs(clinical_fasta, exist_ok=True)
        dirs["clinical_fasta"] = clinical_fasta
    return dirs

# find existing .fna and .gff files in the reference directory
def find_existing_reference_files(reference_dir: str) -> tuple:
    existing_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))
    existing_gff = glob.glob(os.path.join(reference_dir, "*.gff"))
    return existing_fasta, existing_gff
