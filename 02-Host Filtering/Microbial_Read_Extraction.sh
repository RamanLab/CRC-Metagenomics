#!/bin/bash

# CONFIGURATION (Modify these paths as needed)

INDEX="path/to/hg38_index" 
INPUT_DIR="path/to/trimmed_fastq"
OUTPUT_DIR="path/to/non_human_reads"
LOG_DIR="path/to/logs"

# Create output and log directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# PROCESSING LOOP

# Safely loop through R1 files
for r1_file in "${INPUT_DIR}"/trimmed_*_R1.fastq.gz; do
    
    # Check if any files actually matched the glob pattern
    [ -e "$r1_file" ] || continue

    # Extract Sample ID from filename (e.g., trimmed_SampleA_R1.fastq.gz -> SampleA)
    filename=$(basename "$r1_file")
    ID=$(echo "$filename" | sed 's/trimmed_//;s/_R1.fastq.gz//')
        
    # Define corresponding R2 path
    R2="${INPUT_DIR}/trimmed_${ID}_R2.fastq.gz"

    if [ -f "$R2" ]; then
        echo "Processing Sample: $ID"
            
        # Run Bowtie2 mapping
        # Saves alignment summary to the logs folder and filters out non-human reads
        bowtie2 -x "$INDEX" \
            -1 "$r1_file" \
            -2 "$R2" \
            --no-discordant \
            --un-conc-gz "${OUTPUT_DIR}/${ID}_non_human_R%.fastq.gz" \
            > /dev/null \
            2> "${LOG_DIR}/${ID}_bowtie2.log"
            
        echo "Finished processing $ID. Summary saved to log."
    else
        echo "ERROR: Matching R2 file not found for ID $ID"
        echo "Looked for: $R2"
    fi
done

echo "All samples processed successfully."