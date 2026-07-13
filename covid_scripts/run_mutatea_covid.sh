#!/bin/bash

#SBATCH --job-name=mutatea_sarscov2
#SBATCH --output=mutatea_sarscov2_%j.out
#SBATCH --error=mutatea_sarscov2_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Activate conda environment
source ~/.bashrc
conda activate flu_cli

# Change to the mutatea directory
cd /data/tisza/analyses/crm/mutatea

# Run mutatea command
python mutatea.py -p sars_cov2 \
    -m /data/tisza/analyses/crm/mutatea/wastewater_metadata \
    -pr /data/contract/TEPHI/ \
    -ref /data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_covid \
    -c /data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_covid \
    -o /data/tisza/analyses/crm/mutatea/covid_scripts \
    -f -l -a

echo "mutatea covid job completed"
