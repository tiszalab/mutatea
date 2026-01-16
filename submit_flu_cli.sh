#!/bin/bash
#SBATCH --job-name=covid_flu_cli
#SBATCH --output=covid_flu_cli_%j.out
#SBATCH --error=covid_flu_cli_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --mem=128G

# exclusive means it will use a whole node for my job, each node has 32 cpus 
# set a cap limit at 128G

# Initialize conda for bash
source /mmfs1/apps/miniconda3/etc/profile.d/conda.sh

# Activate conda environment
conda activate flu_cli

# Ensure varmint is in PATH (installed at ~/.local/bin/varmint)
export PATH="$HOME/.local/bin:$PATH"

# Navigate to flu_cli directory
cd /data/tisza/analyses/crm/flu_cli

# Run flu_cli with COVID data
python flu_cli.py \
  -s Sars-Cov2 \
  -m /data/tisza/analyses/crm/flu_cli/wastewater_metadata \
  -sr /data/tisza/analyses/crm/cli_outdated/covid_cli/clinical_input_data \
  -ref /data/tisza/analyses/crm/cli_outdated/covid_cli/clinical_input_data \
  -c /data/tisza/analyses/crm/cli_outdated/covid_cli/clinical_input_data \
  -o /data/tisza/analyses/crm/flu_cli

# specified the output so it doesn't overwrite the existing run