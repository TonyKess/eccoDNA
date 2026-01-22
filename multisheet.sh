#!/bin/bash
# Define your target directories
for DIR in /path/to/data1 /path/to/data2 /path/to/data3; do
    echo "Processing $DIR..."
    bash create_samplesheet.sh "$DIR"
done