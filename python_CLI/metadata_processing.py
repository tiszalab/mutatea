#!/usr/bin/env python

###################### SETUP ######################
# load modules
import argparse
import pandas as pd
import glob
import sys, os
import logging
import gzip
import shutil
from Bio import SeqIO
from pathlib import Path
import time

# load in functions from metadata_funcs
# crm: make sure to update names of functions being imported as they change in metadata_funcs.py
#try:
#    from .metadata_funcs import load_metadata, add_region, add_month_year, unzip_reference_files, split_clinical_fasta
#except:
#    from metadata_funcs import load_metadata, add_region, add_month_year, unzip_reference_files, split_clinical_fasta

# convert string to boolean for argparse
def str2bool(x):
    if isinstance(x, bool):
       return x
    if x.lower() in ("yes", "y"):
        return True
    elif x.lower() in ('no', 'n'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# start timer
start_time = time.perf_counter()

# let the user know what this pipeline does and the files required
print("\nThis is a pipeline for processing different serotypes of Influenza A in wastewater data, choosing a reference genome, and (optionally) pulling and processing matched clinical data. \nRequired inputs include: \n - wastewater metadata \n - reference FASTA \n - reference GFF\nThis pipeline will walk you through getting these files and processing them\n")

# request flu subtype
subtype = input("Enter the flu subtype you want to analyze (H1N1, H3N2, or H5N1): ").strip() 

# reprompt the user if it's not one of the set flu subtypes, make sure all subtypes are uppercase for NCBI search
while subtype.upper() not in ["H1N1", "H3N2", "H5N1"]:
    subtype = input("Please enter one of the following flu subtypes (H1N1, H3N2, or H5N1): ").strip()

# set default file path for output to be the folder where this python script is kept
output_path_default = os.path.dirname(os.path.abspath(__file__))

# ask user where they want their output directory (where the processed input files and outputted alignment files will be saved)
default_output_dir = os.path.join(output_path_default, f"{subtype.upper()}_align")

output_dir = input(
    "\nCreate a directory to save the processed input files and the outputted alignment files\n"
    f"(If you hit enter, you can choose the default: {default_output_dir}): "
).strip(" '\"")

# created a default output path for myself
if output_dir == "":
    output_dir = default_output_dir

# make sure the directory exists
os.makedirs(output_dir, exist_ok=True)



###################### WASTEWATER METADATA PROCESSING ######################
# create a subfolder in the output directory for the cleaned metadata files
metadata_dir = os.path.join(output_dir, "metadata_files")

# create the output directory if it doesn't exist
if not os.path.exists(metadata_dir):
    os.makedirs(metadata_dir)

# ask if the user also wants to process clinical data
run_clinical = input("\nWill you also want to process Influenza A reads from Texas clinical data of the same time range as the inputted wastewater data? (y/n): ").strip()
while run_clinical.upper() not in ["Y", "YES", "", "N", "NO"]:
    run_clinical = input("Please enter either y or n: ").strip()

# added in default for myself
if run_clinical == "":
    run_clinical = "y"

# request the file path of the metadata folder
metadata_folder = input("\nEnter the file path of your folder containing the metadata xlsx files: ").strip(" '\"")

# added in default for myself
if metadata_folder == "":
    metadata_folder = (f"{output_path_default}/wastewater_metadata")

# test to see if the metadata folder exists
while not os.path.exists(metadata_folder):
    metadata_folder = input("Enter the file path of your folder containing the metadata xlsx files:").strip(" '\"")

# pull out the xlsx files from the folder path given
metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))

# load in metadata files
md_list=[pd.read_excel(file) for file in metadata_files]

# merge metadata files into one large metadata file
# got idea from https://stackoverflow.com/questions/20908018/import-multiple-excel-files-into-python-pandas-and-concatenate-them-into-one-dat
metadata=pd.concat(md_list, ignore_index=True)

# ask user if they want their wastewater data labelled by public health region
region_request = input("\nDo you also want your wastewater data labelled by public health region? This can be useful for later visualization of how mutations are spreading (y/n): ").strip() 
while region_request.upper() not in ["Y", "YES", "", "N", "NO"]:
    subtype = input("Please enter either y or n: ").strip()

# added in default for myself
if region_request == "":
    region_request = "Y"

# create a dictionary of expected cities and their public health regions
if region_request.upper() in ["Y", "YES"]:
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

    # use the dictionary to add a "Region" column to themetadata
    metadata["Region"] = metadata["City"].map(city_region)

    # add a warning for any unexpected cities, that way the user will know if they need to update the city_region dictionary
    if len(metadata.loc[metadata["Region"].isna(), "City"].unique()) > 0:
        print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())
    else:
        print("All cities in the metadata were successfully assigned to public health regions!")
else:
    pass

# add a column for month_year to the metadata
metadata["Month_Year"] = metadata["Date"].dt.strftime("%m.%Y")

# added this so I can send them test metadata instead of them running the whole dataset
# issue is that the earlier metadata does not include a SiteCode column
if "SiteCode" not in metadata.columns:
    metadata["SiteCode"] = pd.NA

# organization of the metadata depends on the user input
if region_request.upper() in ["Y", "YES"]:
    metadata = metadata[
        [
            "City",
            "Sample_ID",
            "Site",
            "Date",
            "Flow",
            "PoolID",
            "SiteCode",
            "Region",
            "Month_Year",
        ]
    ]
else:
    metadata = metadata[
        [
            "City",
            "Sample_ID",
            "Site",
            "Date",
            "Flow",
            "PoolID",
            "SiteCode",
            "Month_Year",
        ]
    ]

# export metadata as a tsv in the output directory 
metadata.to_csv(f"{metadata_dir}/metadata_wastewater_combined.csv", sep=",", index=False)

# find the time range of the wastewater samples
earliest_date = metadata["Date"].min()
latest_date = metadata["Date"].max()

# print the time range of the wastewater samples
print(f"The wastewater samples range from {earliest_date.strftime('%m/%d/%Y')} to {latest_date.strftime('%m/%d/%Y')}\, you should try to use clinical data that matches this time range")

# I looked at how the dates were formatted in the email and used this website https://strftime.org/
start_str = earliest_date.strftime("%Y-%m-%dT00:00:00.00Z")
end_str = latest_date.strftime("%Y-%m-%dT23:59:59.00Z")

###################### CLINICAL METADATA PROCESSING ######################
# added this if statement so clinical data is only processed if the user said yes
if run_clinical.upper() in ["Y", "YES"]:
    # request file path of clinical metadata
    clinical_metadata_path = input("\nPlease enter the file path of your clinical metadata csv: ").strip(" '\"")

    # add test to confirm they gave the path of a csv
    while (not clinical_metadata_path.endswith(".csv")) or (not os.path.isfile(clinical_metadata_path)):
        clinical_metadata_path = input("Please enter the file path of your clinical metadata csv: ").strip(" '\"")

    # load in clinical metadata
    clinical_metadata = pd.read_csv(clinical_metadata_path)

    # reformat dates in clinical metadata
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")

    # export metadata as tsv
    clinical_metadata.to_csv(f"{metadata_dir}/metadata_clinical_{subtype}.tsv", sep="\t", index=False)


# moved clinical fasta input line to improve run time
    # request file path of clinical fasta
    clinical_fasta_path = input("After downloading the clinical fasta, please enter the file path of your clinical fasta: ").strip(" '\"")

    # add test to confirm they gave the path of a fasta
    while (not clinical_fasta_path.endswith(".fasta")) or (not os.path.isfile(clinical_fasta_path)):
        print("Error: please enter a valid file path that ends in .fasta")
        clinical_fasta_path = input("Enter the file path of your clinical fasta: ").strip(" '\"")



###################### CHOOSE RELEVANT LINEAR REFERENCE GENOME ######################
# offer the same reference genomes I've been using
default_ref = input("\nDo you want to use the default reference genome? (y/n): ").strip()
while default_ref.upper() not in ["Y", "YES", "", "N", "NO"]:
    default_ref = input("Please enter either y or n: ").strip()

# added in default for myself
if default_ref == "":
    default_ref = "Y"

# create a subfolder in the output directory for the cleaned reference files (if it doesn't already exist)
reference_dir = os.path.join(output_dir, "reference_files")

# give the NCBI Virus link for the default reference genome of respective flu subtype if they don't have the files
if not os.path.exists(reference_dir) and default_ref.upper() in ["Y", "YES"]:
    if subtype.lower() == "h1n1":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/308/895/GCA_039308895.1_ASM3930889v1/")
    elif subtype.lower() == "h3n2":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/301/835/GCA_039301835.1_ASM3930183v1/")
    elif subtype.lower() == "h5n1":
        print("Note that H5N1 is not common in the United States, this is the reference genome I found that was closest: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/465/435/GCA_039465435.1_ASM3946543v1/")

###################### PROCESS FILES OF REFERENCE GENOME ######################
# create the output directory if it doesn't exist
if not os.path.exists(reference_dir):
    os.makedirs(reference_dir)

# try to auto-load existing reference files if the user wants the default reference
if default_ref.upper() in ["Y", "YES"] and os.path.isdir(reference_dir):
    # crm: I only load in unzipped versions because zipped files would not be in the output directory
    existing_fasta = glob.glob(os.path.join(reference_dir, "*.fna"))
    existing_gff = glob.glob(os.path.join(reference_dir, "*.gff"))

    # if the user wants default reference files and they were used in a previous run, use them
    if existing_fasta and existing_gff:
        path_ref_fasta = existing_fasta[0].strip(" '\"")
        path_ref_gff = existing_gff[0].strip(" '\"")
        print(f"\nUsing existing reference FASTA: {path_ref_fasta}")
        print(f"\nUsing existing reference GFF: {path_ref_gff}")

    # if only the fasta is downloaded
    elif existing_fasta and not existing_gff:
        path_ref_fasta = existing_fasta[0].strip(" '\"")
        if subtype.lower() == "h1n1":
            print("No gff was found in your folder, download the genomic.gff.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/308/895/GCA_039308895.1_ASM3930889v1/")
        elif subtype.lower() == "h3n2":
            print("No gff was found in your folder, download the genomic.gff.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/301/835/GCA_039301835.1_ASM3930183v1/")
        elif subtype.lower() == "h5n1":
            print("No gff was found in your folder, download the genomic.gff.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/465/435/GCA_039465435.1_ASM3946543v1/")
        path_ref_gff = input("After downloading the reference gff, please enter the file path of your gff.gz or gff: ").strip(" '\"")

    # if only the gff is downloaded
    elif existing_gff and not existing_fasta:
        path_ref_gff = existing_gff[0].strip(" '\"")
        if subtype.lower() == "h1n1":
            print("No fna was found in your folder, download the genomic.fna.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/308/895/GCA_039308895.1_ASM3930889v1/")
        elif subtype.lower() == "h3n2":
            print("No fna was found in your folder, download the genomic.fna.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/301/835/GCA_039301835.1_ASM3930183v1/")
        elif subtype.lower() == "h5n1":
            print("No fna was found in your folder, download the genomic.fna.gz from here: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/465/435/GCA_039465435.1_ASM3946543v1/")
        path_ref_fasta = input("After downloading the reference fasta, please enter the file path of your fna.gz or fna: ").strip(" '\"")
        
# load in the reference fasta
if default_ref.upper() in ["N", "NO"]:
    path_ref_fasta = input("After downloading the reference fasta, please enter the file path of your fna.gz or fna: ").strip(" '\"")

# added the .fna option in case the user unzips the file themselves
while (not path_ref_fasta.endswith(".fna.gz")) and (not path_ref_fasta.endswith(".fna")) or (not os.path.isfile(path_ref_fasta)):
    print("Error: please enter a valid file path that ends in .fna.gz or .fna")
    path_ref_fasta = input("Enter the file path of your reference fasta, make sure the file name ends in fna.gz or fna: ").strip(" '\"")

# unzip if the file path ends in fna.gz
if path_ref_fasta.endswith(".fna.gz"):

    # set the output path as the reference_files folder in the output directory
    fasta_name = os.path.basename(path_ref_fasta[:-3])
    ref_fasta = os.path.join(reference_dir, fasta_name)

    # unzip the fna.gz file to the output directory
    with gzip.open(path_ref_fasta, "rb") as f_in, open(ref_fasta, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"\nUnzipped FASTA to: {ref_fasta}\n")
else:
    # the file is already unzipped if the file name ends in .fna
    fasta_name = os.path.basename(path_ref_fasta)
    ref_fasta = os.path.join(reference_dir, fasta_name)
    if path_ref_fasta != ref_fasta:
        shutil.copy(path_ref_fasta, ref_fasta)
    
# load in reference gff
if default_ref.upper() in ["Y", "YES"] and not os.path.exists(reference_dir):
    path_ref_gff = input("After downloading the reference gff, please enter the file path of your reference gff.gz or gff: ").strip(" '\"")

# added the .gff option in case the user unzips the file themselves
while (not path_ref_gff.endswith(".gff.gz")) and (not path_ref_gff.endswith(".gff")) or (not os.path.isfile(path_ref_gff)):
    print("Error: please enter a valid file path that ends in .gff.gz or .gff")
    path_ref_gff = input("\nEnter the file path of your reference gff, make sure the file name ends in gff.gz or gff: ").strip(" '\"")

# unzip if the file path ends in gff.gz
if path_ref_gff.endswith(".gff.gz"):

    # set the output path to the reference_files folder in the output directory
    gff_name = os.path.basename(path_ref_gff[:-3])
    ref_gff = os.path.join(reference_dir, gff_name)

    # unzip the gff.gz file to the output directory
    with gzip.open(path_ref_gff, "rb") as f_in, open(ref_gff, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"\nUnzipped GFF to: {ref_gff}\n")

else:
    # the file is already unzipped if the file name ends in .gff
    gff_name = os.path.basename(path_ref_gff)
    ref_gff = os.path.join(reference_dir, gff_name)
    if path_ref_gff != ref_gff:
        shutil.copy(path_ref_gff, ref_gff)



###################### SPLIT CLINICAL FASTA BY MONTH ######################
# creating monthly clinical lists to split the clinical fasta later
if run_clinical.upper() in ["Y", "YES"]:
    # create a subfolder in the output directory for the monthly clinical lists that will be created to sort the clinical fasta later
    clinical_output = os.path.join(output_dir, "clinical_output")

    # create the output directory if it doesn't exist
    if not os.path.exists(clinical_output):
        os.makedirs(clinical_output)

    # create a subfolder in the output directory for the monthly clinical lists that will be created to sort the clinical fasta later
    clinical_lists = os.path.join(clinical_output, "monthly_lists")

    # create the output directory if it doesn't exist
    if not os.path.exists(clinical_lists):
        os.makedirs(clinical_lists)
    
    # creates a list of accessions for each month
    for month, group in clinical_metadata.groupby("Month_Year"):
        out_path = Path(clinical_lists) / f"{month}_list.txt"
        group["Accession"].to_csv(out_path, index=False, header=False)
    print(f"\nSorted clinical {subtype} accessions by month")

    # create a subfolder in the output directory for the monthly clinical fastas that will be created later
    clinical_fasta = os.path.join(clinical_output, "monthly_fasta")

    # create the output directory if it doesn't exist
    if not os.path.exists(clinical_fasta):
        os.makedirs(clinical_fasta)

# load in the clinical fasta and split by the monthly lists
if run_clinical.upper() in ["Y", "YES"]:
    # idea from https://stackoverflow.com/questions/51990373/loop-through-dictionary-to-match-dictionary-key-to-a-list-of-values-and-append-d
    # load clinical fasta as a dictionary with the accessiona as the key
    records_by_id = SeqIO.to_dict(SeqIO.parse(clinical_fasta_path, "fasta"))

    # loop through the monthly lists and split the clinical fasta
    # I got the glob idea from https://docs.python.org/3/library/glob.html
    # for every txt file in my clinical_lists directory
    for list_file in Path(clinical_lists).glob("*.txt"):
        # read in the list
        with open(list_file) as f:
            accessions = f.read().splitlines()
        
        # get all accessions for each month
        month_accessions = [records_by_id[a] for a in accessions if a in records_by_id]

        # get the month_year from the file name of the list
        month_year = list_file.name.split("_")[0]
        
        # create the output directory if it doesn't exist, there was an issue when I had the os.path.exists check here
        clinical_monthly_fasta = f"{clinical_fasta}/{month_year}.fasta"

        # export clinical fasta by month, idea from https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        SeqIO.write(month_accessions, clinical_monthly_fasta, "fasta")

    print(f"\nThe clinical FASTA file of {subtype} has been split by month")



###################### FIND AND PULL WASTEWATER READS ######################









###################### STATS FOR CRM ######################
# end timer
end_time = time.perf_counter()

# If I ran it (wastewater metadata was loaded from my default path), tell me how long my script took
if metadata_folder == "/mmfs1/home/u255582/CascadeProjects/flu_mutatome_pipelines/python_CLI/wastewater_metadata":
    print(f"\nTime taken: {end_time - start_time:.2f} seconds\n")

# crm: keep getting error here saying "metadata_processing is not defined"
# crm: call the main function
#if __name__ == "__main__":
#    metadata_processing()