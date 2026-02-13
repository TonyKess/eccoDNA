#!/bin/bash

# Define the target directory for containers, defaults to where ecoDNAflow is 
CONTAINER_DIR="./containers"
mkdir -p $CONTAINER_DIR

echo "Downloading Apptainer images to $CONTAINER_DIR..."

# Array of image names and their URLs
# Format: "filename|url"
images=(
    "pear.0.9.6.3.sif|https://depot.galaxyproject.org/singularity/pear:0.9.6--3"
    "vsearch.2.21.1.sif|https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0"
    "seqkit.2.3.0.sif|https://depot.galaxyproject.org/singularity/seqkit:2.3.0--h9ee0642_0"
    "blast.2.1.3.sif|https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0"
    "pandas.1.5.2.sif|https://depot.galaxyproject.org/singularity/pandas:1.5.2"
)

for entry in "${images[@]}"; do
    IFS="|" read -r FILENAME URL <<< "$entry"
    
    if [ ! -f "$CONTAINER_DIR/$FILENAME" ]; then
        echo "Pulling $FILENAME..."
        apptainer pull "$CONTAINER_DIR/$FILENAME" "$URL"
    else
        echo "$FILENAME already exists, skipping."
    fi
done

echo "Setup complete. Images are located in $CONTAINER_DIR"
