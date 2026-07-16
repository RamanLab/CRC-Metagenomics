#!/bin/bash

# --- Configuration ---
DB="PathtoKraken2NTDatabase"
IN_DIR="Pathto/NonHumanReads"
OUT_DIR="Pathto/KrakenResultsNTDatabasePoint9"
THREADS=32


# Create directory structure
mkdir -p "$OUT_DIR/reports" "$OUT_DIR/bracken" "$OUT_DIR/logs"


# --- 1. Classification & Re-estimation Loop ---
for R1 in "$IN_DIR"/*_R1.fastq.gz; do
    # Identify the corresponding R2 file
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    
    echo "Processing Sample: $SAMPLE" | tee -a "$OUT_DIR/logs/pipeline.log"
    
    # Run Kraken2
    kraken2 --db "$DB" \
            --paired "$R1" "$R2" \
            --threads "$THREADS" \
            --confidence 0.9 \
            --use-names \
            --report "$OUT_DIR/reports/${SAMPLE}.report" \
            --output - >> "$OUT_DIR/logs/pipeline.log" 2>&1

    # Run Bracken
    bracken -d "$DB" \
            -i "$OUT_DIR/reports/${SAMPLE}.report" \
            -o "$OUT_DIR/bracken/${SAMPLE}.bracken" \
            -l S >> "$OUT_DIR/logs/pipeline.log" 2>&1
done

# --- 2. Obtain Merged TSV Results ---
echo "Merging all sample results into a single abundance table..." | tee -a "$OUT_DIR/logs/pipeline.log"

# This uses the Bracken utility script to create a combined matrix
combine_bracken_outputs.py \
    --files "$OUT_DIR/bracken"/*.bracken \
    --output "$OUT_DIR/NII_nt_point9_Abundance.tsv"

echo "Pipeline Complete at $(date)" | tee -a "$OUT_DIR/logs/pipeline.log"
echo "Final abundance table: $OUT_DIR/NII_nt_point9_Abundance.tsv is obtained"
