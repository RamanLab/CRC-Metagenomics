#!/bin/bash

# --- Configuration ---
DB="PathtoKraken2ntDB"
IN_DIR="PathTo/AssemblyFastaFiles"
OUT_DIR="PathTo/KrakenResultsAssemblyNTDatabasePoint9"
THREADS=32

# Create directory structure
mkdir -p "$OUT_DIR/reports" "$OUT_DIR/logs"

# --- Classification Loop ---
for SAMPLE_DIR in "$IN_DIR"/*/; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    FASTA="$SAMPLE_DIR/final_assembly.fasta"

    if [[ ! -f "$FASTA" ]]; then
        echo "Skipping $SAMPLE — no final_assembly.fasta found" | tee -a "$OUT_DIR/logs/pipeline.log"
        continue
    fi

    echo "Processing Sample: $SAMPLE" | tee -a "$OUT_DIR/logs/pipeline.log"

    kraken2 --db "$DB" \
            "$FASTA" \
            --threads "$THREADS" \
            --confidence 0.9 \
            --use-names \
            --report "$OUT_DIR/reports/${SAMPLE}.report" \
            --output - >> "$OUT_DIR/logs/pipeline.log" 2>&1
done


MERGED="$OUT_DIR/NII_assembly_kraken_merged.tsv"
echo -e "sample\tpercent\tclade_reads\ttaxon_reads\trank\ttaxid\tname" > "$MERGED"

for REPORT in "$OUT_DIR/reports/"*.report; do
    SAMPLE=$(basename "$REPORT" .report)
    awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"} {print sample, $1, $2, $3, $4, $5, $6}' "$REPORT" >> "$MERGED"
done

echo "Pipeline Complete at $(date)" | tee -a "$OUT_DIR/logs/pipeline.log"
echo "Final merged table: $MERGED"
