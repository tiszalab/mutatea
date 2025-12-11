#!/usr/bin/env python

###################### SETUP ######################
# load modules
import pandas as pd
import gzip
import shutil
from pathlib import Path
from Bio import SeqIO

###################### SETUP ######################
# load in and merge metadata files
def load_metadata(metadata_folder):

# add public health region column to merged metadata
def add_region(metadata):

# add month_year column to metadata
def add_month_year(metadata):

# unzip reference files if needed
def unzip_reference_files(reference_fasta, reference_gff):

# split clinical fasta by month
def split_clinical_fasta(clinical_fasta, output_dir):
