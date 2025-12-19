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

###### crm: this is maybe(?) functional
# process reference files
def process_reference_file(input_folder: str, reference_dir: str) -> list:
    # make sure reference_dir exists
    os.makedirs(reference_dir, exist_ok=True)
    output_paths=[]

    # find all .fna.gz and .gff.gz files in the input folder and unzip into the reference_dir
    gz_files = (
        glob.glob(os.path.join(input_folder, "*fna.gz"))
        + glob.glob(os.path.join(input_folder, "*gff.gz"))
    )
    for gz_file in gz_files:
        # Get the output filename by removing .gz extension
        filename = os.path.basename(gz_file)[:-3]
        out_path = os.path.join(reference_dir, filename)
        with gzip.open(gz_file, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"\nUnzipped to: {out_path}\n")
        output_paths.append(out_path)

    # copy any uncompressed reference files into the reference_dir
    for pattern in ("*.fna", "*.gff"):
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
    existing_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))
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
        clinical_fasta_month = os.path.join(output_dir, f"{month_year}.fasta")
        SeqIO.write(month_accessions, clinical_fasta_month, "fasta")

# crm: confirm correct setup of the include_region input here
# create output directories
def create_output_directories(output_dir: str, subtype: str, include_region:bool = True, include_clinical: bool = False) -> dict:
    dirs = {}
    # main output directory with subtype-specific subfolder
    dirs["output"] = os.path.join(output_dir, f"{subtype}_align")
    os.makedirs(dirs["output"], exist_ok=True)

    # cleaned and merged metadata directory
    dirs["metadata_dir"] = os.path.join(dirs["output"], "metadata_files")
    os.makedirs(dirs["metadata_dir"], exist_ok=True)
    
    # output file directory
    dirs["reference_dir"] = os.path.join(dirs["output"], "reference_files")
    os.makedirs(dirs["reference_dir"], exist_ok=True)

    # create file directory for alignment files
    dirs["alignment_dir"] = os.path.join(dirs["output"], "alignment_files")
    os.makedirs(dirs["alignment_dir"], exist_ok=True)

    # create subfolders in the alignment directory
    # crm: we could also just not create the ww folder if no clinical analysis is done?
    # crm confirm this is valid structure for include_clinical
    dirs["wastewater_dir"] = os.path.join(dirs["alignment_dir"], "wastewater")
    os.makedirs(dirs["wastewater_dir"], exist_ok=True)
    
    # wastewater lists
    dirs["wastewater_lists_dir"] = os.path.join(dirs["wastewater_dir"], "lists")
    os.makedirs(dirs["wastewater_lists_dir"], exist_ok=True)

    # different structure depending on if region is included
    if include_region:
        # create subfolder for lists_month_region
        dirs["wastewater_list_reg"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month_reg")
        os.makedirs(dirs["wastewater_list_reg"], exist_ok=True)

        # create subfolder for lists_month
        dirs["wastewater_list_month"] = os.path.join(dirs["wastewater_lists_dir"], "lists_month")
        os.makedirs(dirs["wastewater_list_month"], exist_ok=True)
    else:
        pass
    # crm: check that pass works here, I only want to create those subfolders if region was included
        
    
    # wastewater merged bams
    dirs["merged_bams"] = os.path.join(dirs["wastewater_dir"], "merged_bams")
    os.makedirs(dirs["merged_bams"], exist_ok=True)

    # create merged_bams subfolders
    # wastewater bams merged by month
    dirs["merged_bams_month"] = os.path.join(dirs["merged_bams"], "merged_bams_month")        
    os.makedirs(dirs["merged_bams_month"], exist_ok=True)

    # wastewater bams merged by month and region
    if include_region:
        dirs["merged_bams_month_region"] = os.path.join(dirs["merged_bams"], "merged_bams_month_region")
        os.makedirs(dirs["merged_bams_month_region"], exist_ok=True)
        
    # create subfolder for alignment files by pool
    dirs["pools"] = os.path.join(dirs["wastewater_dir"], "pools")
    os.makedirs(dirs["pools"], exist_ok=True)

    # optional: clinical directories
    if include_clinical:
        # parent folder for clinical output
        dirs["clinical_output"] = os.path.join(dirs["alignment_dir"], "clinical_output")
        os.makedirs(dirs["clinical_output"], exist_ok=True)

        # folder for the lists of accessions by month
        dirs["clinical_lists_month"] = os.path.join(dirs["clinical_output"], "clinical_lists_month")
        os.makedirs(dirs["clinical_lists_month"], exist_ok=True)

        # crm: want to later remove these clinical fasta files, could maybe save them to a temp dir?
        # folder for the clinical fastas split by month
        dirs["clinical_fasta_month"] = os.path.join(dirs["clinical_output"], "clinical_fasta_month")
        os.makedirs(dirs["clinical_fasta_month"], exist_ok=True)

        # folder for the clinical bam files that were merged by month
        dirs["clinical_bam_month"] = os.path.join(dirs["clinical_output"], "clinical_bam_month")
        os.makedirs(dirs["clinical_bam_month"], exist_ok=True)       

    # create tsv_output folder to later catch tsv files
    dirs["tsv_output"] = os.path.join(dirs["output"], "tsv_output")
    os.makedirs(dirs["tsv_output"], exist_ok=True)

    # create the directories I described
    return dirs