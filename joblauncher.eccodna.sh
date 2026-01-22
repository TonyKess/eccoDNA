#!/bin/bash

# Loop through every config file in the configs directory
for CONFIG_PATH in configs/project_*.config; do
    echo "Submitting job for $CONFIG_PATH..."
    
    # Use --export to pass the specific path to the sbatch environment
    sbatch --export=ALL,CONFIG_FILE="$CONFIG_PATH" submit.ednanf.customconfig.sh
done