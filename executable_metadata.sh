#!/bin/bash

# This script ensures the conda environment is activated before running the Python script

# Get conda base path and source it
source $(conda info --base)/etc/profile.d/conda.sh

# Activate the flu_cli environment
conda activate flu_cli

# Check if activation was successful
if [ $? -eq 0 ]; then
    echo "The conda environment 'flu_cli' was successfully activated"
    echo "Running metadata_processing.py..."
    echo ""
    # Change to the script directory and run the Python script
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    cd "$SCRIPT_DIR"
    python metadata_processing.py
else
    echo "Failed to activate conda environment 'flu_cli'"
    exit 1
fi
