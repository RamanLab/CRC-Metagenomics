#!/bin/bash
# Activate Environment
source PathToConda/bin/activate
conda activate humanN


# Path Definitions
INPUT_DIR="/PathTo/NonHumanReads"
OUTPUT_DIR="/PathTo/HumanNResults"
DB_DIR="/PathTo/HumannDB"

#Get Sample Name (Using SGE_TASK_ID for all samples)
# This reads the line number corresponding to the current job task
SAMPLE_NAME=$(sed -n "${SGE_TASK_ID}p" /PathTo/SampleList.txt)

# Join R1 and R2
cat ${INPUT_DIR}/${SAMPLE_NAME}_non_human_R1.fastq.gz ${INPUT_DIR}/${SAMPLE_NAME}_non_human_R2.fastq.gz > /PathTo/JoinedFiles/${SAMPLE_NAME}_joined.fastq.gz

# Run HUMAnN
# Using 10 threads as per your sample command
humann --input /PathTo/JoinedFiles/${SAMPLE_NAME}_joined.fastq.gz \
  --output ${OUTPUT_DIR}/${SAMPLE_NAME} \
  --threads 10 \
  --remove-temp-output \
  --metaphlan-options "--bowtie2db ${DB_DIR}/metaphlan_db --index mpa_vJun23_CHOCOPhlAnSGB_202403"

# Cleanup the temporary joined file
rm /PathTo/JoinedFiles/${SAMPLE_NAME}_joined.fastq.gz
