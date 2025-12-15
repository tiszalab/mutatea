import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
from subprocess import Popen
import argparse
import glob
import os

# load in and merge metadata files
def load_metadata(metadata_folder:str):
    try:
        metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))
        md_list=[pd.read_excel(file) for file in metadata_files]
        return metadata=pd.concat(md_list, ignore_index=True)
    except:
        

# add month_year column to metadata
def add_month_year(metadata:str):
    metadata["Month_Year"] = metadata["Date"].dt.strftime("%m.%Y")

# optional: add public health region column to merged metadata
def add_region(metadata:str):
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
    # use the dictionary to add a "Region" column to the metadata
    metadata["Region"] = metadata["City"].map(city_region)

    # add a warning for any unexpected cities, that way the user will know if they need to update the city_region dictionary
    if len(metadata.loc[metadata["Region"].isna(), "City"].unique()) > 0:
        print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())
    else:
        print("All cities in the metadata were successfully assigned to public health regions!")


# crm: feel like I could use the same function for both of these, need to check to remove redundancy
# unzip reference files if needed
def unzip_reference_fasta(path_ref_fasta:str, unzipref:str):
    return Popen(['gunzip', '-c', path_ref_fasta, '>', unzipref], stdout=(unzipref, "w"), stderr=STDOUT)

def unzip_reference_gff(path_ref_gff:str, unzipref:str):
    return Popen(['gunzip', '-c', path_ref_gff, '>', unzipref], stdout=(unzipref, "w"), stderr=STDOUT)

# split clinical fasta by month
def split_clinical_fasta(clinical_fasta:str, output_dir:str):

