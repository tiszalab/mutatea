import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
import glob
import os

# crm: need filter to only merge files that have specific columns
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

# get the time range from the metadata
def get_date_range(metadata: pd.DataFrame) -> tuple:
    earliest_date = metadata["Date"].min()
    latest_date = metadata["Date"].max()
    return earliest_date, latest_date

# export the metadata as a csv to output path
def export_metadata(metadata: pd.DataFrame, metadata_dir: str, sep: str = ",") -> None:
    metadata.to_csv(f"{metadata_dir}/metadata_wastewater_combined.csv", sep=sep, index=False)

# load in clinical metadata and add column for month_year
def load_clinical_metadata(clinical_metadata_path: str) -> pd.DataFrame:
    clinical_metadata = pd.read_csv(clinical_metadata_path)
    # make collection date a datetime object
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    # add month_year column to the clinical metadata
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")
    return clinical_metadata

# load in clinical fasta
def load_clinical_fasta(clinical_fasta_path: str) -> list:
    clinical_fasta = SeqIO.parse(clinical_fasta_path, "fasta")
    return clinical_fasta

###### crm: this is maybe(?) functional
# process reference files
def process_reference_file(input_folder: str, reference_dir: str) -> list:
    # make sure reference_dir exists
    os.makedirs(reference_dir, exist_ok=True)
    output_paths=[]

    # find all .gz files in the input folder and unzip into the reference_dir
    gz_files = glob.glob(os.path.join(input_folder, "*.gz"))
    for gz_file in gz_files:
        # Get the output filename by removing .gz extension
        filename = os.path.basename(gz_file)[:-3]
        out_path = os.path.join(reference_dir, filename)
        with gzip.open(gz_file, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"\nUnzipped to: {out_path}\n")
        output_paths.append(out_path)

    # copy any uncompressed reference files into the reference_dir
    for pattern in ("*.fasta", "*.fa", "*.fna", "*.gff"):
        for src_path in glob.glob(os.path.join(input_folder, pattern)):
            filename = os.path.basename(src_path)
            out_path = os.path.join(reference_dir, filename)
            # make sure the reference file is not already in the reference_dir
            if os.path.exists(out_path):
                print(f"Reference file {filename} already exists in {reference_dir}.\n")
                continue
            if src_path != out_path:
                shutil.copy(src_path, out_path)
            print(f"Copied to: {out_path}\n")
            output_paths.append(out_path)
    return output_paths

# find existing .fna and .gff files in the reference directory
def find_existing_reference_files(reference_dir: str) -> tuple:
    existing_fasta = (
        glob.glob(os.path.join(reference_dir, "*.fna"))
        + glob.glob(os.path.join(reference_dir, "*.fasta"))
        + glob.glob(os.path.join(reference_dir, "*.fa"))
    )
    existing_gff = glob.glob(os.path.join(reference_dir, "*.gff"))
    return existing_fasta, existing_gff

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

    # optional: clinical directories
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