###################### SETUP ######################
# Confirm the user has all required modules installed
print("Please check the module_requirements.txt file and install any missing modules with pip install -r module_requirements.txt")

# load modules
import pandas as pd
import glob
import os
import gzip
import shutil
from Bio import SeqIO
from pathlib import Path

# let the user know what this pipeline does and the files required
print("\nThis is a pipeline for processing different serotypes of Influenza A in wastewater data, choosing a reference genome, and (optionally) pulling and processing matched clinical data. \nRequired inputs include: \n - wastewater metadata \n - reference FASTA \n - reference GFF\nThis pipeline will walk you through getting these files and processing them\n")

# request flu subtype
subtype = input("Enter the flu subtype you want to analyze (H1N1, H3N2, or H5N1): ").strip() 

# reprompt the user if it's not one of the set flu subtypes
while subtype not in ["H1N1", "H3N2", "H5N1"]:
    subtype = input("Please enter one of the following flu subtypes (H1N1, H3N2, or H5N1): ").strip()

# I got the idea to work with os.path from https://www.reddit.com/r/learnpython/comments/11xybcz/how_can_i_find_the_path_to_the_downloads_folder/
# following os.path functions were chose with the information from https://docs.python.org/3/library/os.path.html
# set default output_path
output_path_default = os.path.expanduser("~/Downloads/python_CLI")

# ask user where they want their output directory (where the processed input files and outputted alignment files will be saved)
output_path = input(
    "\nCreate a directory to save the processed input files and the outputted alignment files\n"
    f"(If you hit enter, you can choose the default: {output_path_default}/{subtype}_align): ").strip()

# I had an issue with copying the filepath from MacOS where sometimes it would include a space or single quote
# I got this idea from chatgpt: removes spaces, single quotes, and double quotes if they exist in the file path
output_path = output_path.strip(" '\"")

# created a default output path for myself
if output_path == "":
    output_path = f"{output_path_default}/{subtype}_align"
    output_path = output_path.strip(" '\"")

# create an output directory called {subtype}_align
output_dir = os.path.join(output_path)

# create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



###################### WASTEWATER METADATA PROCESSING ######################
# create a subfolder in the output directory for the cleaned metadata files
metadata_dir = os.path.join(output_dir, "metadata_files")

# create the output directory if it doesn't exist
if not os.path.exists(metadata_dir):
    os.makedirs(metadata_dir)

# ask if the user also wants to process clinical data
run_clinical = input("\nWill you also want to process Influenza A reads from Texas clinical data of the same time range as the inputted wastewater data? (y/n): ").strip()
while run_clinical not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
    run_clinical = input("Please enter either y or n: ").strip()

# added in default for myself
if run_clinical == "":
    run_clinical = "y"

# request the file path of the metadata folder
metadata_folder = input("\nEnter the file path of your folder containing the metadata xlsx files: ").strip()

# remove spaces, single quotes, and double quotes if they exist in the file path
metadata_folder = metadata_folder.strip(" '\"")

# test to see if the metadata folder exists
while not os.path.exists(metadata_folder):
    metadata_folder = input("Enter the file path of your folder containing the metadata xlsx files:").strip()
    metadata_folder = metadata_folder.strip(" '\"")

# pull out the xlsx files from the folder path given
metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))

# load in metadata files
md_list=[pd.read_excel(file) for file in metadata_files]

# merge metadata files into one large metadata file
# got idea from https://stackoverflow.com/questions/20908018/import-multiple-excel-files-into-python-pandas-and-concatenate-them-into-one-dat
metadata=pd.concat(md_list, ignore_index=True)

# ask user if they want their wastewater data labelled by public health region
region_request = input("\nDo you also want your wastewater data labelled by public health region? This can be useful for later visualization of how mutations are spreading (y/n): ").strip() 
while region_request not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
    subtype = input("Please enter either y or n: ").strip()

# added in default for myself
if region_request == "":
    region_request = "y"

# create a dictionary of expected cities and their public health regions
if region_request.lower() in ["y", "Y", "yes", "Yes"]:
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
if region_request.lower() in ["y", "Y", "yes", "Yes"]:
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

# print the time range of the wastewater samples
earliest_date = metadata["Date"].min()
latest_date = metadata["Date"].max()
time_match = input("\nDo you want to know the time range of the wastewater samples? (y/n): ").strip() 
while time_match not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
    time_match = input("Please enter either y or n: ").strip()

# added in default for myself
if time_match == "":
    time_match = "y"

# print the time range of the wastewater samples if the user requests it
if time_match.lower() in ["y", "Y", "yes", "Yes"]:
    print(f"The wastewater samples range from {earliest_date.strftime("%m/%d/%Y")} to {latest_date.strftime("%m/%d/%Y")}")
    if run_clinical.lower() in ["y", "yes"]:
        print("You should try to use clinical data that matches this time range")

# reformat dates for NCBI Virus URL (specifically need to remove the spaces so the url works)
# I looked at how the dates were formatted in the email and used this website https://strftime.org/
start_str = earliest_date.strftime("%Y-%m-%dT00:00:00.00Z")
end_str = latest_date.strftime("%Y-%m-%dT23:59:59.00Z")

# offer NCBI link and instructions for getting clinical files if the user doesn't already have them downloaded
if run_clinical.lower() in ["y", "yes"]:
    have_clinical_files = input("\nDo you already have the clinical metadata (csv) and clinical reads (FASTA) downloaded? (y/n): ").strip()
    while have_clinical_files not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
        have_clinical_files = input("Please enter either y or n: ").strip()
    
    # added in default for myself
    if have_clinical_files == "":
        have_clinical_files = "y"

    # ask if the user wants an NCBI link if they don't have the clinical files
    if have_clinical_files.lower() in ["n", "no"]:
        want_ncbi_help = input("\nDo you want an NCBI Virus link and instructions for downloading the clinical metadata and FASTA? (y/n): ").strip()
        while want_ncbi_help not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
            want_ncbi_help = input("Please enter either y or n: ").strip()
        
        # added in default for myself
        if want_ncbi_help == "":
            want_ncbi_help = "n"

        # offer NCBI link and instructions for getting required clinical files if the user doesn't already have them
        if want_ncbi_help.lower() in ["y", "yes"]:
            print(f"\nHere is the link: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&GenomeCompleteness_s=complete&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&CollectionDate_dr={start_str}%20TO%20{end_str}&Serotype_s={subtype}&USAState_s=TX")
            print(f"\nYou will need to:\n - Download all records as a nucleotide FASTA \n - Download the metadata as a csv and select all, making sure to include the accession with version \n")



###################### CLINICAL METADATA PROCESSING ######################
# added this if statement so clinical data is only processed if the user said yes
if run_clinical.lower() in ["y", "yes"]:
    # request file path of clinical metadata
    clinical_metadata_path = input("\nPlease enter the file path of your clinical metadata csv: ").strip()

    # remove spaces, single quotes, and double quotes if they exist in the file path
    clinical_metadata_path = clinical_metadata_path.strip(" '\"")

    # add test to confirm they gave the path of a csv
    while (not clinical_metadata_path.endswith(".csv")) or (not os.path.isfile(clinical_metadata_path)):
        clinical_metadata_path = input("Please enter the file path of your clinical metadata csv: ").strip()
        clinical_metadata_path = clinical_metadata_path.strip(" '\"")

    # load in clinical metadata
    clinical_metadata = pd.read_csv(clinical_metadata_path)

    # reformat dates in clinical metadata
    # idea from https://pandas.pydata.org/docs/reference/api/pandas.to_datetime.html
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")

    # export metadata as tsv
    clinical_metadata.to_csv(f"{metadata_dir}/metadata_clinical_{subtype}.tsv", sep="\t", index=False)



###################### CHOOSE RELEVANT LINEAR REFERENCE GENOME ######################
# offer the same reference genomes I've been using
default_ref = input("\nDo you want to use the default reference genome? (y/n): ").strip()
while default_ref not in ["y", "Y", "yes", "Yes", "", "n", "N", "no", "No"]:
    default_ref = input("Please enter either y or n: ").strip()

# added in default for myself
    if default_ref == "":
        default_ref = "y"

# give the NCBI Virus link for the default reference genome of respective flu subtype
if default_ref.lower() in ["y", "Y", "yes", "Yes"]:
    # chatgpt helped me troubleshoot this part (apparently the subtype needed to be lowercase)
    if subtype.lower() == "h1n1":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/308/895/GCA_039308895.1_ASM3930889v1/")
    elif subtype.lower() == "h3n2":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/301/835/GCA_039301835.1_ASM3930183v1/")
    elif subtype.lower() == "h5n1":
        print("Note that H5N1 is not common in the United States, this is the reference genome I found that was closest: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/465/435/GCA_039465435.1_ASM3946543v1/")

# offer the option to use their own reference genome and give them the NCBI Virus link
else:
    print(f"\nHere is the link: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Genome&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&GenomeCompleteness_s=complete&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&CollectionDate_dr={start_str}%20TO%20{end_str}&Serotype_s={subtype}&USAState_s=TX")
    print(f"You will want to pick a reference genome from around the beginning of your time range, which is {start_str} \nYou will need to download the GFF and FASTA of the selected genome")

# tell the user which required files they need to download for the reference genome
print("\nDownload the files for genomic.gff.gz and genomic.fna.gz")



###################### PROCESS FILES OF REFERENCE GENOME ######################
# create a subfolder in the output directory for the cleaned reference files
reference_dir = os.path.join(output_dir, "reference_files")

# create the output directory if it doesn't exist
if not os.path.exists(reference_dir):
    os.makedirs(reference_dir)

# load in the reference fasta
path_ref_fasta = input("After downloading the reference fasta, please enter the file path of your fna.gz or fna: ").strip()

# remove spaces, single quotes, and double quotes if they exist in the file path
path_ref_fasta = path_ref_fasta.strip(" '\"")

# added the .fna option in case the user unzips the file themselves
# chatgpt helped me fix this issue, I needed to put an "and" between the file type options
while (not path_ref_fasta.endswith(".fna.gz")) and (not path_ref_fasta.endswith(".fna")) or (not os.path.isfile(path_ref_fasta)):
    print("Error: please enter a valid existing file path that ends in .fna.gz or .fna")
    path_ref_fasta = input("Enter the file path of your reference fasta, make sure the file name ends in fna.gz or fna: ").strip()
    path_ref_fasta = path_ref_fasta.strip(" '\"")

# unzip if the file path ends in fna.gz
if path_ref_fasta.endswith(".fna.gz"):

    # set the output path as the reference_files folder in the output directory
    fasta_name = os.path.basename(path_ref_fasta[:-3])
    ref_fasta = os.path.join(reference_dir, fasta_name)

    # chatgpt recommended using shutil to copy the content to the new file
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
    print(f"\nFASTA file is already unzipped: {ref_fasta}\n")

# load in reference gff
path_ref_gff = input("After downloading the reference gff, please enter the file path of your reference gff.gz or gff: ").strip()

# remove spaces, single quotes, and double quotes if they exist in the file path
path_ref_gff = path_ref_gff.strip(" '\"")

# added the .gff option in case the user unzips the file themselves
while (not path_ref_gff.endswith(".gff.gz")) and (not path_ref_fasta.endswith(".gff")) or (not os.path.isfile(path_ref_gff)):
    print("Error: please enter a valid existing file path that ends in .gff.gz or .gff")
    path_ref_gff = input("\nEnter the file path of your reference gff, make sure the file name ends in gff.gz or gff: ").strip()
    path_ref_gff = path_ref_gff.strip(" '\"")

# unzip if the file path ends in gff.gz
if path_ref_gff.endswith(".gff.gz"):

    # set the output path to the reference_files folder in the output directory
    gff_name = os.path.basename(path_ref_gff[:-3])
    ref_gff = os.path.join(reference_dir, gff_name)

    # chatgpt recommended using shutil to copy the content to the new file
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
    print(f"\nGFF file is already unzipped: {ref_gff}\n")



###################### SPLIT CLINICAL FASTA BY MONTH ######################
# creating monthly clinical lists to split the clinical fasta later
if run_clinical.lower() in ["y", "yes"]:
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
    
    # chatgpt helped me with this loop, creates a list of accessions for each month
    for month, group in clinical_metadata.groupby("Month_Year"):
        out_path = Path(clinical_lists) / f"{month}_list.txt"
        group["Accession"].to_csv(out_path, index=False, header=False)
    print(f"Sorted clinical {subtype} accessions by month")

    # create a subfolder in the output directory for the monthly clinical fastas that will be created later
    clinical_fasta = os.path.join(clinical_output, "monthly_fasta")

    # create the output directory if it doesn't exist
    if not os.path.exists(clinical_fasta):
        os.makedirs(clinical_fasta)

# load in the clinical fasta and split by the monthly lists
if run_clinical.lower() in ["y", "yes"]:
    # request file path of clinical fasta
    clinical_fasta_path = input("After downloading the clinical fasta, please enter the file path of your clinical fasta: ").strip()

    # remove spaces, single quotes, and double quotes if they exist in the file path
    clinical_fasta_path = clinical_fasta_path.strip(" '\"")

    # add test to confirm they gave the path of a fasta
    while (not clinical_fasta_path.endswith(".fasta")) or (not os.path.isfile(clinical_fasta_path)):
        print("Error: please enter a valid existing file path that ends in .fasta")
        clinical_fasta_path = input("Enter the file path of your clinical fasta: ").strip()
        clinical_fasta_path = clinical_fasta_path.strip(" '\"")

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
        
        # get all accessions for each month, chatgpt helped me with this line 
        month_accessions = [records_by_id[a] for a in accessions if a in records_by_id]

        # get the month_year from the file name of the list
        month_year = list_file.name.split("_")[0]
        
        # create the output directory if it doesn't exist, there was an issue when I had the os.path.exists check here
        clinical_monthly_fasta = f"{clinical_fasta}/{month_year}.fasta"

        # export clinical fasta by month, idea from https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        SeqIO.write(month_accessions, clinical_monthly_fasta, "fasta")

    print(f"\nThe clinical FASTA file of {subtype} has been split by month\n")