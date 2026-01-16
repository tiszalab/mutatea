#!/bin/bash
#SBATCH --job-name=test_varmint
#SBATCH --output=test_varmint.out

export PATH="$HOME/.local/bin:$PATH"
which varmint
varmint --help