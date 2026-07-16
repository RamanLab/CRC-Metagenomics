#!/bin/bash
cd PathTo/NonHumanReads

echo $$
conda init bash
. PathToConda/bin/activate
PATH=PathToMetaWRAP/bin:$PATH

conda activate metawrap-env


metawrap assembly         -1 PathTo/NonHumanReads/SampleID_R1.fastq.gz         -2 PathTo/NonHumanReads/SampleID_R2.fastq.gz         -m 200         -t 40         -o PathTo/SampleID
