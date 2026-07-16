#!/bin/bash

# Create directories for outputs (relative to where you run the script)
mkdir -p trimmed_fastq
mkdir -p fastp_reports

echo "Starting batch adapter removal..."

# Loop through all R1 files using the full path
for r1_file in PathToYourFolder/*_R1*.fastq.gz
do
    # To get the filename without the path for output naming
    r1_name=$(basename "$r1_file")
    
    # Identify the matching R2 file path
    # This replaces _R1 with _R2 in the full path
    r2_file=${r1_file/_R1/_R2}
    r2_name=$(basename "$r2_file")

    # Checks if R2 exists
    if [ ! -f "$r2_file" ]; then
        echo "Skip: Could not find matching R2 for $r1_name"
        continue
    fi

    # Creates a clean sample name 
    sample_name=$(echo "$r1_name" | cut -d'_' -f1)

    echo "---------------------------------------------------"
    echo "Processing Sample: $sample_name"

    # Run fastp
    ./fastp \
        -i "$r1_file" \
        -I "$r2_file" \
        -o "trimmed_fastq/trimmed_${r1_name}" \
        -O "trimmed_fastq/trimmed_${r2_name}" \
        --detect_adapter_for_pe \
        --cut_tail \
        -h "fastp_reports/${sample_name}.html" \
        -j "fastp_reports/${sample_name}.json" \
        --thread 12

done

echo "==================================================="
echo "Batch trimming complete!"
echo "Clean files are in: $(pwd)/trimmed_fastq/"
