# Confirm the user has all required modules installed
print("Please check the module_requirements.txt file and install any missing modules with pip install -r module_requirements.txt")

# load modules
import pandas as pd
import glob
import os
import gzip
import shutil
import subprocess


# let the user know what this pipeline does and the files required
### crm: adjust the output line depending on what I get done
## crm: clean this intro line, it's messy
print("\nThis is a pipeline for taking reads of different serotypes of Influenza A in Texas wastewater data and aligning them to a reference genome to output a tsv of the nonsynonymous mutations occurring in each source. \nRequired inputs include: \n - wastewater metadata \n - wastewater FASTA files \n - reference FASTA \n - reference GFF\nThis pipeline will walk you through getting these files and processing them\n")

# request flu subtype
subtype = input("Enter the flu subtype you want to analyze (H1N1, H3N2, or H5N1): ").strip() 

# add in test to reprompt the user if it's not one of the set flu subtypes
while subtype not in ["H1N1", "H3N2", "H5N1"]:
    subtype = input("Please enter one of the following flu subtypes (H1N1, H3N2, or H5N1): ").strip()


# ask user where they want their output directory (where the processed input files, outputted alignment files, and tsv of mutations will be saved )
# default output_path
output_path_default = f"/Users/camillemazurek2025/python_CLI"

output_path = input(
    "\nCreate a directory to save the processed input files, the outputted alignment files, and the tsv of mutations\n"
    f"(The default will be: {output_path_default}/{subtype}_align):"
).strip()


## crm: make sure you can explain the difference between output_path and output_dir here, looks repetitive
if output_path == "":
    output_path = f"{output_path_default}/{subtype}_align"

# remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
output_path = output_path.strip(" '\"")

# create output directory called {subtype}_align
output_dir = os.path.join(output_path)

# create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



########## WASTEWATER METADATA PROCESSING ##########
# create a subfolder in the output directory for the cleaned metadata files
metadata_dir = os.path.join(output_dir, "metadata_files")

# create the output directory if it doesn't exist
if not os.path.exists(metadata_dir):
    os.makedirs(metadata_dir)


# ask if the user also wants to process clinical data
run_clinical = input("\nWill you also want to process Influenza A reads from Texas clinical data of the same time range as the inputted wastewater data? (y/n): ").strip()
while run_clinical not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
    run_clinical = input("Please enter either y or n: ").strip()

# load in metadata
# request file path of metadata
## crm: need to add the input for metadata folder later
# metadata_folder = input("Enter the file path of your metadata xlsx files: ").strip()
# add test to see if this is a valid file path
## crm: need to set to run before submitting
# while metadata_folder not in os.listdir():
#    metadata_folder = input("Enter the file path of your metadata xlsx files: ").strip()
## remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
#metadata_folder = metadata_folder.strip(" '\"")

metadata_folder="/Users/camillemazurek2025/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/data2/metadata"
metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))

# load in metadata files
md_list=[pd.read_excel(file) for file in metadata_files]

# merge metadata files into one large metadata file
### crm: not sure if I should ignore the index here
metadata=pd.concat(md_list, ignore_index=True)

# create a dictionary of expected cities and their public health regions
# ask user if they want their wastewater data split by public health region
region_request = input("\nDo you also want your wastewater data split by public health region? This can be useful for later visualization of how mutations are spreading (y/n): ").strip() 
while region_request not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
    subtype = input("Please enter either y or n: ").strip()
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

    # use dictionary to add region column to metadata
    metadata["Region"] = metadata["City"].map(city_region)

    # add warning for any unexpected cities, would need to update city_region dictionary
    if len(metadata.loc[metadata["Region"].isna(), "City"].unique()) > 0:
        print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())
    else:
        print("All cities in the metadata were successfully assigned to public health regions!")
else:
    metadata["Region"] = "All"


# add a column for month_year to the metadata
metadata["Month_Year"] = metadata["Date"].dt.strftime("%m.%Y")

# organize the metadata in a specific way
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

# export metadata as tsv
metadata.to_csv(f"{metadata_dir}/metadata_wastewater_combined.csv", sep=",", index=False)


# print the time range of the wastewater samples
earliest_date = metadata["Date"].min()
latest_date = metadata["Date"].max()
time_match = input("\nDo you want to know the time range of the wastewater samples? (y/n): ").strip() 
while time_match not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
    time_match = input("Please enter either y or n: ").strip()
if time_match.lower()  in ["y", "Y", "yes", "Yes"]:
    print(f"The wastewater samples range from {earliest_date.strftime("%m/%d/%Y")} to {latest_date.strftime("%m/%d/%Y")}, you should use clinical data that matches this time range")

# reformat dates for NCBI Virus URL (specifically need to remove the spaces so the url works)
start_str = earliest_date.strftime("%Y-%m-%dT00:00:00.00Z")
end_str = latest_date.strftime("%Y-%m-%dT23:59:59.00Z")

# offer NCBI link and instructions for getting clinical files if the user doesn't already have them downloaded
if run_clinical.lower() in ["y", "yes"]:
    have_clinical_files = input("\nDo you already have the clinical metadata (csv) and clinical reads (FASTA) downloaded? (y/n): ").strip()
    while have_clinical_files not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
        have_clinical_files = input("Please enter either y or n: ").strip()

    if have_clinical_files.lower() in ["n", "no"]:
        want_ncbi_help = input("\nDo you want an NCBI link and instructions for downloading the clinical metadata and FASTA? (y/n): ").strip()
        while want_ncbi_help not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
            want_ncbi_help = input("Please enter either y or n: ").strip()

        if want_ncbi_help.lower() in ["y", "yes"]:
            print(f"\nHere is the link: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&GenomeCompleteness_s=complete&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&CollectionDate_dr={start_str}%20TO%20{end_str}&Serotype_s={subtype}&USAState_s=TX")
            print(f"\nYou will need to:\n - Download all records as a nucleotide FASTA \n - Download the metadata as a csv and select all, making sure to include the accession with version \n")




########## CLINICAL METADATA PROCESSING ##########
# added this if statement so clinical data is only processed if the user said yes
if run_clinical.lower() in ["y", "yes"]:
    # request file path of clinical metadata
    clinical_metadata_path = input("\nPlease enter the file path of your clinical metadata csv: ").strip()

    # remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
    clinical_metadata_path = clinical_metadata_path.strip(" '\"")

    # add test to confirm they gave the path of a csv
    while (not clinical_metadata_path.endswith(".csv")) or (not os.path.isfile(clinical_metadata_path)):
        print("Error: please enter a valid existing file path that ends in .csv")
        clinical_metadata_path = input("Please enter the file path of your clinical metadata csv: ").strip()


    # load in clinical metadata
    clinical_metadata = pd.read_csv(clinical_metadata_path)

    # reformat dates in clinical metadata
    clinical_metadata["Collection_Date"] = pd.to_datetime(clinical_metadata["Collection_Date"], errors="coerce")
    clinical_metadata["Month_Year"] = clinical_metadata["Collection_Date"].dt.strftime("%m.%Y")

    # export metadata as tsv
    clinical_metadata.to_csv(f"{metadata_dir}/metadata_clinical_{subtype}.tsv", sep="\t", index=False)


########## CHOOSE RELEVANT LINEAR REFERENCE GENOME ##########

# offer the same reference genomes I've been using
default_ref = input("\nDo you want to use the default reference genome? (y/n): ").strip()
while default_ref not in ["y", "Y", "yes", "Yes", "n", "N", "no", "No"]:
    default_ref = input("Please enter either y or n: ").strip()
if default_ref.lower() in ["y", "Y", "yes", "Yes"]:
    # crm: chatgpt helped me troubleshoot this part (subtype needs to be lowercase)
    if subtype.lower() == "h1n1":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/308/895/GCA_039308895.1_ASM3930889v1/")
    elif subtype.lower() == "h3n2":
        print("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/301/835/GCA_039301835.1_ASM3930183v1/")
    elif subtype.lower() == "h5n1":
        print("Note that H5N1 is not common in the United States, this is the reference genome I found that was closest: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/465/435/GCA_039465435.1_ASM3946543v1/")

# offer the option to use their own reference genome
#### crm: maybe let them see the genomes +- 6 months of the earliest date in the time range
else:
    print(f"\nHere is the link: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Genome&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&GenomeCompleteness_s=complete&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&CollectionDate_dr={start_str}%20TO%20{end_str}&Serotype_s={subtype}&USAState_s=TX")
    print(f"You will want to pick a reference genome from around the beginning of your time range, which is {start_str} \n You will need to download the GFF and FASTA of the selected genome")

print("\nDownload the files for genomic.gff.gz and genomic.fna.gz")







########## PROCESS FILES OF REFERENCE GENOME ##########
# create a subfolder in the output directory for the cleaned reference files
reference_dir = os.path.join(output_dir, "reference_files")

# create the output directory if it doesn't exist
if not os.path.exists(reference_dir):
    os.makedirs(reference_dir)

# load in reference fasta
path_ref_fasta = input("After downloading the reference fasta, please enter the file path of your fna.gz or fna: ").strip()

# remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
path_ref_fasta = path_ref_fasta.strip(" '\"")

# added the .fna option in case the user unzips the file themselves
# crm: chatgpt helped me realize I needed to put an "and" between the file type options
while (not path_ref_fasta.endswith(".fna.gz")) and (not path_ref_fasta.endswith(".fna")) or (not os.path.isfile(path_ref_fasta)):
    print("Error: please enter a valid existing file path that ends in .fna.gz or .fna")
    path_ref_fasta = input("Enter the file path of your reference fasta, make sure the file name ends in fna.gz or fna: ").strip()

# unzip if the file path ends in fna.gz
if path_ref_fasta.endswith(".fna.gz"):

    # set the output path as the same place as the input path but without the .gz extension
    fasta_name = os.path.basename(path_ref_fasta[:-3])
    ref_fasta = os.path.join(reference_dir, fasta_name)

    ## crm: chatgpt recommended using shutil
    # unzip the fna.gz file to a new FASTA file in the same place as the input file
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

# remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
path_ref_gff = path_ref_gff.strip(" '\"")

# added the .gff option in case the user unzips the file themselves
while (not path_ref_gff.endswith(".gff.gz")) and (not path_ref_fasta.endswith(".gff")) or (not os.path.isfile(path_ref_gff)):
    print("Error: please enter a valid existing file path that ends in .gff.gz or .gff")
    path_ref_gff = input("Enter the file path of your reference gff, make sure the file name ends in gff.gz or gff: ").strip()

# unzip if the file path ends in gff.gz
if path_ref_gff.endswith(".gff.gz"):
    

    ### crm: want to correct this to download to output_dir
    # set the output path as the same place as the input path but without the .gz extension
    gff_name = os.path.basename(path_ref_gff[:-3])
    ref_gff = os.path.join(reference_dir, gff_name)

    ## crm: chatgpt recommended using shutil
    # unzip the gff.gz file to a new gff file in the same place as the input file
    with gzip.open(path_ref_gff, "rb") as f_in, open(ref_gff, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"\nUnzipped GFF to: {ref_gff}\n")
else:
    # the file is already unzipped if the file name ends in .gff
    gff_name = os.path.basename(path_ref_gff)
    ref_gff = os.path.join(reference_dir, gff_name)
    ## crm: know the mechanism for the shutil copy line
    if path_ref_gff != ref_gff:
        shutil.copy(path_ref_gff, ref_gff)
    print(f"\nGFF file is already unzipped: {ref_gff}\n")


########## SPLIT CLINICAL FASTA BY MONTH ##########
# need to include if statement so clinical data is only processed if the user said yes
if run_clinical.lower() in ["y", "yes"]:
    # request file path of clinical fasta
    clinical_fasta_path = input("After downloading the clinical fasta, please enter the file path of your clinical fasta: ").strip()

    # remove quotes if they exist in the file path (was an issue when copying file path on MacOS)
    clinical_fasta_path = clinical_fasta_path.strip(" '")

    # add test to confirm they gave the path of a fasta
    while (not clinical_fasta_path.endswith(".fasta")) or (not os.path.isfile(clinical_fasta_path)):
        print("Error: please enter a valid existing file path that ends in .fasta")
        clinical_fasta_path = input("Enter the file path of your clinical fasta: ").strip()

    # load in clinical fasta
    clinical_fasta = pd.read_csv(clinical_fasta_path, header=None, names=["Sequence"])







    # split clinical fasta by month
    #clinical_fasta_by_month = clinical_fasta.groupby("Month_Year")

    # export clinical fasta by month
    ## crm: set to run before submitting
    # clinical_fasta_by_month_filepath = input("Enter the file path where you want your processed clinical fasta by month saved: ").strip() 
    # while clinical_fasta_by_month_filepath not in os.listdir():
    #    clinical_fasta_by_month_filepath = input("Enter the file path where you want your processed clinical fasta by month saved: ").strip()
    # clinical_fasta_by_month.to_csv(f"{clinical_fasta_by_month_filepath}/{subtype}_clinical_fasta_by_month.tsv", sep="\t", index=False)

    #print(f"\nExporting the processed clinical fasta by month as {subtype}_clinical_fasta_by_month.tsv\n")
    #clinical_fasta_by_month.to_csv(f"/Users/camillemazurek2025/Downloads/{subtype}_clinical_fasta_by_month.tsv", sep="\t", index=False)





########## LOAD IN WASTEWATER READS ##########
print("\nNow we will align the wastewater reads to the reference genome and merge BAMs by month.\nIf you chose to also split the wastewater data by public health region, that will also be done here.")





# ask user for one or more PoolIDs to process
pool_ids = input(
    "Enter one or more PoolIDs to process (comma-separated, matching the metadata PoolID column): "
).strip().split(",")

shell_script = "/Users/camillemazurek2025/flu_mutatome_pipelines/09.04_mergedbam_novarm.sh"

for raw_pool in pool_ids:
    poolid = raw_pool.strip()
    if not poolid:
        continue
    print(f"\nRunning alignment and BAM merging for subtype {subtype} and pool {poolid}...\n")
    try:
        subprocess.run(
            ["bash", shell_script, subtype, poolid],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running script for pool {poolid}: {e}")