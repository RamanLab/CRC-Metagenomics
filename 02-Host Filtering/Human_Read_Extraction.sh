#!/bin/bash

# ENVIRONMENT SETUP

source /PathTo/conda/etc/profile.d/conda.sh
conda activate bowtie2

# Verify tools are available
for tool in bowtie2 samtools picard bedtools; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool not found. Check conda environment."
        exit 1
    fi
done

echo ">>> Environment is ready to work with"


# PATHS & VARIABLES

# Indexes
HG38_INDEX="PathTohg38Index"
BAC_INDEX="PathToBacterialIndex"

# Input
FASTQ_DIR="PathToTrimmedFastq"

# Output
OUT_DIR="PathToOutputDirectory"

# Reference files
BED_OF_NONPAR="PathToNonPARBedFile"

# Threads
THREADS= NumberOfThreadsToUse 


# If running as SGE task array, get ID from list
SAMPLE_LIST="/data/durga/NII/sample_list.txt"

if [ ! -z "$SGE_TASK_ID" ]; then
    ID=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_LIST")
    echo ">>> SGE Task $SGE_TASK_ID → Sample: $ID"
fi

# Validate ID is set
if [ -z "$ID" ]; then
    echo "ERROR: ID not set. Either SGE_TASK_ID or -v ID= required."
    exit 1
fi


# Construct input paths
R1="${FASTQ_DIR}/trimmed_${ID}_R1.fastq.gz"
R2="${FASTQ_DIR}/trimmed_${ID}_R2.fastq.gz"

# Validate input FASTQs exist
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: Input FASTQs not found for $ID"
    echo "  Expected R1: $R1"
    echo "  Expected R2: $R2"
    exit 1
fi

# Output directories
SAMPLE_DIR="${OUT_DIR}/${ID}"
BAM_DIR="${SAMPLE_DIR}/BAM"
METRIC_DIR="${SAMPLE_DIR}/metrics"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${BAM_DIR}" "${METRIC_DIR}" "${LOG_DIR}"

echo ">>> Starting sample: $ID at $(date)"
echo ">>> Using $THREADS threads"


# Check if final BAM already exists — skip if so
FINAL_BAM="${BAM_DIR}/non_bacterial_${ID}_mapped.bam"

if [ -f "$FINAL_BAM" ]; then
    echo ">>> Final BAM already exists for $ID. Skipping..."
    exit 0
fi

# COUNT ORIGINAL READS (NUM0)
# Count read pairs in original trimmed FASTQ
echo ">>> Counting input reads for $ID..."

NUM0=$(zcat "$R1" | wc -l | awk '{print $1/4}')

echo ">>> NUM0 (trimmed input pairs): $NUM0"



# ALIGN TO HUMAN GENOME (hg38)
# Output: raw SAM with ALL reads (mapped + unmapped)

echo ">>> Aligning $ID to hg38 at $(date)"

bowtie2 \
    -p $THREADS \
    --no-discordant \
    -x ${HG38_INDEX} \
    -1 "$R1" \
    -2 "$R2" \
    -S "${BAM_DIR}/${ID}.sam" \
    2> "${METRIC_DIR}/bowtie2_hg38_${ID}.log"

# Check bowtie2 succeeded
if [ $? -ne 0 ]; then
    echo "ERROR: Bowtie2 hg38 alignment failed for $ID"
    exit 1
fi

echo ">>> hg38 alignment done at $(date)"


# CONVERT, SORT, AND FILTER TO PROPERLY PAIRED MAPPED
# 
# -f 3 breakdown:
#   bit 1 (0x1) = read is paired
#   bit 2 (0x2) = read is properly paired (both mapped, correct orientation)
#   Combined = keep ONLY properly paired mapped reads
#
# What gets REMOVED here:
#   - Unmapped reads 
#   - Mixed pairs (one maps, one doesn't)
#   - Discordant pairs
#   - Multimappers in wrong orientation

echo ">>> Converting, sorting, filtering to properly paired at $(date)"

samtools view -Sb -@ $THREADS "${BAM_DIR}/${ID}.sam" | \
samtools sort -@ $THREADS -m 3G | \
samtools view -@ $THREADS -u -f 3 \
    -o "${BAM_DIR}/${ID}_mapped_pre.bam"

# Check pipeline succeeded
if [ $? -ne 0 ]; then
    echo "ERROR: SAM to filtered BAM conversion failed for $ID"
    exit 1
fi

# Remove SAM immediately — large file
rm "${BAM_DIR}/${ID}.sam"

samtools index "${BAM_DIR}/${ID}_mapped_pre.bam"

# Count properly paired mapped read PAIRS
NUM1=$(samtools view -c "${BAM_DIR}/${ID}_mapped_pre.bam" | awk '{print $1/2}')

echo ">>> NUM1 (properly paired mapped pairs after hg38): $NUM1"
echo ">>> Filtering done at $(date)"



# ADD READ GROUPS
# Picard MarkDuplicates requires read group tags (RG)
# Without these it either fails or processes incorrectly
#
# RGID uses actual sample ID (better than hardcoded FLOWCELLID as in reference 

echo ">>> Adding read groups for $ID at $(date)"

picard -Xmx100g AddOrReplaceReadGroups \
    I="${BAM_DIR}/${ID}_mapped_pre.bam" \
    O="${BAM_DIR}/${ID}_mapped.bam" \
    RGID="${ID}" \
    RGLB="${ID}_library_1" \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM="${ID}" \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

if [ $? -ne 0 ]; then
    echo "ERROR: AddOrReplaceReadGroups failed for $ID"
    exit 1
fi

# Remove pre-RG BAM
rm "${BAM_DIR}/${ID}_mapped_pre.bam" \
   "${BAM_DIR}/${ID}_mapped_pre.bam.bai"

echo ">>> Read groups added at $(date)"


# REMOVE PCR DUPLICATES
# Picard MarkDuplicates identifies reads that are PCR duplicates — identical copies from amplification not from independent DNA molecules
# REMOVE_DUPLICATES=true  : physically removes them (not just flags)
# ASSUME_SORTED=true      : skips re-sorting, saves time
#

echo ">>> Removing duplicates for $ID at $(date)"

picard -Xmx100g MarkDuplicates \
    I="${BAM_DIR}/${ID}_mapped.bam" \
    O="${BAM_DIR}/rm_dup_${ID}_mapped.bam" \
    M="${METRIC_DIR}/rm_dup_${ID}_metrics.txt" \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

if [ $? -ne 0 ]; then
    echo "ERROR: MarkDuplicates failed for $ID"
    exit 1
fi

samtools index "${BAM_DIR}/rm_dup_${ID}_mapped.bam"

# Remove intermediate BAMs
rm "${BAM_DIR}/${ID}_mapped.bam" \
   "${BAM_DIR}/${ID}_mapped.bai"

# Count deduplicated read PAIRS (NUM2)
NUM2=$(samtools view -c "${BAM_DIR}/rm_dup_${ID}_mapped.bam" | awk '{print $1/2}')

echo ">>> NUM2 (pairs after duplicate removal): $NUM2"
echo ">>> Duplicate removal done at $(date)"


# BACTERIAL DECONTAMINATION
# Method: Tomofuji et al.
# Screen deduplicated human reads against bacterial database
# Keep only reads that do NOT map to bacteria
#
# Because we applied -f 3 earlier, the rm_dup BAM now contains
# ONLY properly paired human mapped reads going into this step
# This means clean_ids.txt will contain ONLY true human read IDs

echo ">>> Starting bacterial decontamination for $ID at $(date)"

# A. Convert deduplicated human BAM back to FASTQ for screening
picard -Xmx100g SamToFastq \
    I="${BAM_DIR}/rm_dup_${ID}_mapped.bam" \
    F="${BAM_DIR}/rm_dup_${ID}_R1.fastq" \
    F2="${BAM_DIR}/rm_dup_${ID}_R2.fastq" \
    VALIDATION_STRINGENCY=LENIENT

gzip -f "${BAM_DIR}/rm_dup_${ID}_R1.fastq"
gzip -f "${BAM_DIR}/rm_dup_${ID}_R2.fastq"

echo ">>> SamToFastq done at $(date)"

# B. Align against bacterial index
bowtie2 \
    -p $THREADS \
    -x ${BAC_INDEX} \
    -1 "${BAM_DIR}/rm_dup_${ID}_R1.fastq.gz" \
    -2 "${BAM_DIR}/rm_dup_${ID}_R2.fastq.gz" \
    -S "${BAM_DIR}/bac_screen_${ID}.sam" \
    2> "${METRIC_DIR}/bowtie2_bac_${ID}.log"

if [ $? -ne 0 ]; then
    echo "ERROR: Bacterial bowtie2 alignment failed for $ID"
    exit 1
fi

echo ">>> Bacterial alignment done at $(date)"

# C. Extract read IDs where BOTH reads are unmapped to bacteria
#    -f 12  : both R1 AND R2 must be unmapped (flag 4 + flag 8 = 12)
#    -F 256 : exclude secondary alignments
samtools view -@ $THREADS -f 12 -F 256 "${BAM_DIR}/bac_screen_${ID}.sam" | \
    cut -f 1 | sort | uniq > "${BAM_DIR}/clean_ids_${ID}.txt"

CLEAN_ID_COUNT=$(wc -l < "${BAM_DIR}/clean_ids_${ID}.txt")
echo ">>> Clean (non-bacterial) read pairs identified: $CLEAN_ID_COUNT"

# Remove bacterial SAM and temporary FASTQs
rm "${BAM_DIR}/bac_screen_${ID}.sam" \
   "${BAM_DIR}/rm_dup_${ID}_R1.fastq.gz" \
   "${BAM_DIR}/rm_dup_${ID}_R2.fastq.gz"

# D. Filter original human BAM using clean ID list
picard -Xmx100g FilterSamReads \
    I="${BAM_DIR}/rm_dup_${ID}_mapped.bam" \
    O="${FINAL_BAM}" \
    RLF="${BAM_DIR}/clean_ids_${ID}.txt" \
    FILTER=includeReadList \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

if [ $? -ne 0 ]; then
    echo "ERROR: FilterSamReads failed for $ID"
    exit 1
fi

# Remove intermediate files
rm "${BAM_DIR}/rm_dup_${ID}_mapped.bam" \
   "${BAM_DIR}/rm_dup_${ID}_mapped.bam.bai" \
   "${BAM_DIR}/clean_ids_${ID}.txt"

# Count final clean read PAIRS (NUM4)
NUM3=$(samtools view -c "${FINAL_BAM}" | awk '{print $1/2}')

echo ">>> NUM3 (final clean human non-bacterial pairs): $NUM3"
echo ">>> Bacterial decontamination done at $(date)"


# METRICS
echo ">>> Generating metrics for $ID at $(date)"

# 1. Read attrition table
# Tracks how many read pairs survived each step
echo -e "Original_Trimmed\tHuman_Mapped\tDup_Removed\tBac_Removed" \
    > "${METRIC_DIR}/read_attrition_${ID}.txt"
echo -e "${NUM0}\t${NUM1}\t${NUM2}\t${NUM3}" \
    >> "${METRIC_DIR}/read_attrition_${ID}.txt"

# 2. Chromosome-level mapping statistics
samtools idxstats "${FINAL_BAM}" \
    > "${METRIC_DIR}/idxstats_${ID}.txt"

# 4. FastQC on final BAM
conda activate fastqc
mkdir -p "${METRIC_DIR}/fastqc"
fastqc \
    -f bam \
    -t 4 \
    --nogroup \
    "${FINAL_BAM}" \
    -o "${METRIC_DIR}/fastqc/"

echo ">>> All metrics done at $(date)"


# FINAL SUMMARY
echo ""
echo " PIPELINE COMPLETE: $ID"
echo " Final BAM:    ${FINAL_BAM}"
echo " Metrics dir:  ${METRIC_DIR}"
echo ""
echo " Read pairs:"
echo "   Input (trimmed):       $NUM0"
echo "   Human mapped (-f 3):   $NUM1"
echo "   After deduplication:   $NUM2"
echo "   After bac removal:     $NUM3"
echo ""
echo ">>> Pipeline finished for $ID at $(date)"

