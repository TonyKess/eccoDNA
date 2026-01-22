#!/bin/bash
#SBATCH --account=account
#SBATCH --partition=standard
#SBATCH --comment="comment"
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=$PROJECT.ecoDNAflow
#SBATCH --output=logs/edna_%j.out
#SBATCH --error=logs/edna_%j.err
#SBATCH --error=logs/edna_%j.err


# Load Nextflow module (adjust to your cluster)
source ~/.bashrc
mamba activate nf-core

# Set directories for Nextflow and Apptainer
export NXF_TEMP=/tmp
export APPTAINER_TMPDIR=/tmp
export TMPDIR=/tmp

mkdir -p logs

# Define the central pipeline location
PIPELINE_SCRIPT="eDNA_nextflow/"
${CONFIG_FILE:-nextflow.config}

CONFIG_FILE=${CONFIG_FILE:-nextflow.config}

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Configuration file $CONFIG_FILE not found!"
    exit 1
fi

nextflow run main.nf -c "$CONFIG_FILE"

# Run the pipeline
# Nextflow will automatically use the nextflow.config in this directory
nextflow run $PIPELINE_SCRIPT/main.nf \
    -profile slurm \
    -c CONFIG_FILE \
    -resume

echo "Pipeline submitted at $(date)"


