#!/usr/bin/env python

import os

# get all files from input data for e coli
ecoli_dir = "/data/tisza/analyses/crm/mutatea/test_input_data/input_data_ecoli/clinical_samples"
all_files = os.listdir(ecoli_dir)

for filename in all_files:
    # split to get base and read number
    if '_1.fastq' in filename:
        pathogen_name = filename.replace('_1.fastq', '.ecoli_1.fastq')
    elif '_2.fastq' in filename:
        pathogen_name = filename.replace('_2.fastq', '.ecoli_2.fastq')
    else:
        continue

    # rename the file
    old_path = os.path.join(ecoli_dir, filename)
    new_path = os.path.join(ecoli_dir, pathogen_name)
    
    # save new files to ecoli directory
    os.rename(old_path, new_path)