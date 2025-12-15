#!/usr/bin/env python

###################### SETUP ######################
# load modules
import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO
from subprocess import Popen
import argparse
import glob
import os

###################### SETUP ######################
output_path_default = os.path.dirname(os.path.abspath(__file__))
# load in and merge metadata files
def load_metadata(metadata_folder):
    if metadata_folder == "":
        metadata_folder = f"{output_path_default}/wastewater_metadata"
    return glob.glob(os.path.join(metadata_folder,"*.xlsx"))

# add month_year column to metadata
def add_month_year(metadata):

# optional: add public health region column to merged metadata
def add_region(metadata):


# unzip reference files if needed
def unzip_reference_files(path_ref_fasta, path_ref_gff):
    if path_ref_fasta.endswith(".fna.gz"):
        # set the output path as the reference_files folder in the output directory
        fasta_name = os.path.basename(path_ref_fasta[:-3])
        ref_fasta = os.path.join(reference_dir, fasta_name)
        # unzip the fna.gz file to the output directory
        with gzip.open(path_ref_fasta, "rb") as f_in, open(ref_fasta, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    elif path_ref_gff.endswith(".gff.gz"):
        # set the output path as the reference_files folder in the output directory
        gff_name = os.path.basename(path_ref_gff[:-3])
        ref_gff = os.path.join(reference_dir, gff_name)
        with gzip.open(path_ref_gff, "rb") as f_in, open(ref_gff, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        pass
        # crm: need the reference files to be kept in the output directory in the reference_files subfolder
        # not sure if I wrote this script so it only unzips the fasta OR the gff, need to confirm

# split clinical fasta by month
def split_clinical_fasta(clinical_fasta, output_dir):
    return gunzip(clinical_fasta)
