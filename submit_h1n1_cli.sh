#!/bin/bash
#SBATCH --job-name=flu_cli
#SBATCH --output=flu_cli%j.out
#SBATCH --error=flu_cli%j.out
#SBATCH --cpus-per-task=32
#SBATCH --exclusive

# exclusive means it will use a whole node for my job, each node has 32 cpus 

# Initialize conda for bash
source /mmfs1/apps/miniconda3/etc/profile.d/conda.sh

# Activate conda environment
conda activate flu_cli

# Ensure varmint is in PATH (installed at ~/.local/bin/varmint)
export PATH="$HOME/.local/bin:$PATH"

# Navigate to working directory
cd /data/tisza/analyses/crm

# Run flu_CLI with H1N1 data
flu_CLI \
  -s H1N1 \
  -m /data/tisza/analyses/crm/flu_cli/wastewater_metadata \
  -pr /data/service/Pools/EsViritu \
  -ref /data/tisza/analyses/crm/flu_cli/test_input_data/clinical_input_data_H1N1 \
  -c /data/tisza/analyses/crm/flu_cli/test_input_data/clinical_input_data_H1N1 \
  -o /data/tisza/analyses/crm/flu_cli

# specified the output so it doesn't overwrite the existing run