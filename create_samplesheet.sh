#!/bin/bash
# create_samplesheet.sh - Generate sample sheet from FASTQ files

DATA_DIR="$1"  # Directory containing FASTQ files
OUTPUT="samplesheet.csv"

echo "sample,read1,read2" > $OUTPUT

for R1 in ${DATA_DIR}/*_R1*.fastq.gz; do
    R2="${R1/_R1/_R2}"
    SAMPLE=$(basename $R1 | sed 's/_R1.*//')
    
    if [[ -f "$R2" ]]; then
        echo "${SAMPLE},${R1},${R2}" >> $OUTPUT
    else
        echo "Warning: No R2 file found for $R1"
    fi
done

echo "Sample sheet created: $OUTPUT"
cat $OUTPUT
