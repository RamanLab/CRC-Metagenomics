#!/bin/bash

# ENVIRONMENT

source PathToConda/etc/profile.d/conda.sh
conda activate gatk_env



# PATHS & VARIABLES
REF="PathTo/hg38/hg38.fa"
BASE_DIR="PathTo/HUMAN_READ"
OUT_DIR="PathTo/MUTECT2_FINAL_RESULTS"

# Supporting files
GNOMAD="PathTo/GatkEssentials/af-only-gnomad.hg38.vcf.gz"
PON="PathTo/GatkEssentials/1000g_pon.hg38.vcf.gz"
FUNCOTATOR_SOURCES="PathTo/GatkEssentials/funcotator_dataSources.v1.8.hg38.20230908s"

# Threads 
THREADS=21

mkdir -p "${OUT_DIR}/logs"

export REF BASE_DIR OUT_DIR GNOMAD PON FUNCOTATOR_SOURCES THREADS

# All samples
SAMPLES=(30 31 32 34 35 36 38 39 40 41 42 43 45 47 48 49 50 66)


# VALIDATE REQUIRED FILES EXIST BEFORE STARTING
echo ">>> Validating required resource files..."

for f in "$REF" "$GNOMAD" "$PON"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "       Please download all resource files first."
        exit 1
    fi
done

if [ ! -d "$FUNCOTATOR_SOURCES" ]; then
    echo "ERROR: Funcotator sources directory not found: $FUNCOTATOR_SOURCES"
    echo "       Run: gatk FuncotatorDataSourceDownloader --somatic --hg38 --extract-after-download"
    exit 1
fi

echo ">>> All resource files validated"


# HELPER — resolve correct BAM paths per sample --- we used since the naming convention differs for samples 30, 31, 32
get_bam_paths() {
    local i=$1

    # Tumor path — consistent across all samples
    T_DIR="${BASE_DIR}/Tumor-P-${i}/BAM"
    TUMOR_BAM="${T_DIR}/non_bacterial_Tumor-P-${i}_mapped.bam"
    TUMOR_SAMPLE_NAME="Tumor-P-${i}"

    # Non-tumor path — naming differs for 30, 31, 32
    if [[ "$i" == "30" || "$i" == "31" || "$i" == "32" ]]; then
        N_DIR="${BASE_DIR}/Non-tumor-P-${i}/BAM"
        NORMAL_BAM="${N_DIR}/non_bacterial_Non-tumor-P-${i}_mapped.bam"
        NORMAL_SAMPLE_NAME="Non-tumor-P-${i}"
    else
        N_DIR="${BASE_DIR}/Non-Tumor-P-${i}/BAM"
        NORMAL_BAM="${N_DIR}/non_bacterial_Non-Tumor-P-${i}_mapped.bam"
        NORMAL_SAMPLE_NAME="Non-Tumor-P-${i}"
    fi

    export TUMOR_BAM NORMAL_BAM TUMOR_SAMPLE_NAME NORMAL_SAMPLE_NAME
}


# MAIN FUNCTION — full best practices pipeline for one sample
run_bestpractices_pipeline() {
    local i=$1

    echo "============================================"
    echo " Starting Sample P-${i} at $(date)"
    echo "============================================"

    get_bam_paths "$i"

    # Output directories per sample
    SAMPLE_DIR="${OUT_DIR}/Sample_${i}"
    mkdir -p "${SAMPLE_DIR}"

    # Output file paths
    VCF_RAW="${SAMPLE_DIR}/somatic_${i}_raw.vcf.gz"
    VCF_FILTERED="${SAMPLE_DIR}/somatic_${i}_filtered.vcf.gz"
    VCF_FUNCOTATED="${SAMPLE_DIR}/somatic_${i}_funcotated.vcf"
    TSV_OUT="${SAMPLE_DIR}/somatic_${i}_PASS_final.tsv"
    LOG="${SAMPLE_DIR}/log_${i}.txt"

    # Resume check — skip if final TSV already exists
    if [ -f "$TSV_OUT" ]; then
        echo ">>> Sample P-${i} already completed. Skipping..."
        return 0
    fi

    # Validate BAMs exist
    if [ ! -f "$TUMOR_BAM" ]; then
        echo "ERROR: Tumor BAM not found: $TUMOR_BAM"
        return 1
    fi
    if [ ! -f "$NORMAL_BAM" ]; then
        echo "ERROR: Normal BAM not found: $NORMAL_BAM"
        return 1
    fi

    # Auto-index BAMs if needed
    if [ ! -f "${TUMOR_BAM}.bai" ]; then
        echo ">>> Indexing tumor BAM for P-${i}..."
        samtools index -@ $THREADS "$TUMOR_BAM"
    fi
    if [ ! -f "${NORMAL_BAM}.bai" ]; then
        echo ">>> Indexing normal BAM for P-${i}..."
        samtools index -@ $THREADS "$NORMAL_BAM"
    fi

    echo ">>> Tumor BAM:  $TUMOR_BAM"
    echo ">>> Normal BAM: $NORMAL_BAM"


    # Mutect2
    if [ ! -f "$VCF_RAW" ]; then
        echo ">>> STEP 1: Running Mutect2 for P-${i} at $(date)..."

        gatk --java-options "-Xmx160g -XX:ParallelGCThreads=4" Mutect2 \
            -R "$REF" \
            -I "$TUMOR_BAM" \
            -I "$NORMAL_BAM" \
            -tumor "$TUMOR_SAMPLE_NAME" \
            -normal "$NORMAL_SAMPLE_NAME" \
            --germline-resource "$GNOMAD" \
	    --panel-of-normals "$PON" \
            --native-pair-hmm-threads $THREADS \
            -O "$VCF_RAW" \
            2>> "$LOG"

        if [ $? -ne 0 ]; then
            echo "ERROR: Mutect2 failed for P-${i}. Check $LOG"
            return 1
        fi
        echo ">>> STEP 1 done at $(date)"
    else
        echo ">>> STEP 1 already done for P-${i}. Skipping..."
    fi


    # FilterMutectCalls
    if [ ! -f "$VCF_FILTERED" ]; then
        echo ">>> STEP 5: FilterMutectCalls for P-${i} at $(date)..."

        gatk --java-options "-Xmx160g" FilterMutectCalls \
            -R "$REF" \
            -V "$VCF_RAW" \
            -O "$VCF_FILTERED" \
            2>> "$LOG"

        if [ $? -ne 0 ]; then
            echo "ERROR: FilterMutectCalls failed for P-${i}"
            return 1
        fi

        # Count PASS variants
        PASS_COUNT=$(bcftools view -f PASS "$VCF_FILTERED" | grep -v "^#" | wc -l)
        echo ">>> STEP 5 done. PASS variants: $PASS_COUNT at $(date)"
    else
        echo ">>> STEP 5 already done for P-${i}. Skipping..."
    fi


    # Funcotator — Full Biological Annotation

    if [ ! -f "$VCF_FUNCOTATED" ]; then
        echo ">>> STEP 6: Funcotator for P-${i} at $(date)..."

        gatk --java-options "-Xmx160g" Funcotator \
            -R "$REF" \
            -V "$VCF_FILTERED" \
            --ref-version hg38 \
            --data-sources-path "$FUNCOTATOR_SOURCES" \
            --output-file-format VCF \
            -O "$VCF_FUNCOTATED" \
            2>> "$LOG"

        if [ $? -ne 0 ]; then
            echo "ERROR: Funcotator failed for P-${i}"
            return 1
        fi
        echo ">>> STEP 6 done at $(date)"
    else
        echo ">>> STEP 6 already done for P-${i}. Skipping..."
    fi


    # Extract PASS variants to TSV
   
    echo ">>> STEP 7: Extracting PASS variants to TSV for P-${i} at $(date)..."

    echo -e "Sample\tCHROM\tPOS\tREF\tALT\tFILTER\tTumor_AF\tTumor_DP\tTumor_AD\t\
Gene\tVariant_Class\tSecondary_Class\tVariant_Type\t\
cDNA_Change\tProtein_Change\t\
ClinVar_Sig\tClinVar_Disease\tClinVar_RS\t\
COSMIC_Mutations\tCOSMIC_Tissue_Alterations\tCOSMIC_Tissues\t\
dbSNP_RS\tDrugBank\tGO_Bio_Process" \
        > "$TSV_OUT"

    bcftools view -f PASS "$VCF_FUNCOTATED" | \
    bcftools query \
        -s "$TUMOR_SAMPLE_NAME" \
        -f "${i}\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF]\t[%DP]\t[%AD]\t%INFO/FUNCOTATION\n" | \
    awk -F'\t' 'BEGIN{OFS="\t"}{
        funcot = $10
        gsub(/^\[/, "", funcot)
        gsub(/\]$/, "", funcot)
        sub(/,.*/, "", funcot)

        n = split(funcot, arr, "|")

        gene         = (n>=1)   ? arr[1]   : "."
        var_class    = (n>=6)   ? arr[6]   : "."
        sec_class    = (n>=7)   ? arr[7]   : "."
        var_type     = (n>=8)   ? arr[8]   : "."
        cdna_chg     = (n>=17)  ? arr[17]  : "."
        protein_chg  = (n>=19)  ? arr[19]  : "."
        clinvar_sig  = (n>=34)  ? arr[34]  : "."
        clinvar_dis  = (n>=30)  ? arr[30]  : "."
        clinvar_rs   = (n>=44)  ? arr[44]  : "."
        cosmic_mut   = (n>=47)  ? arr[47]  : "."
        cosmic_tot   = (n>=50)  ? arr[50]  : "."
        cosmic_tiss  = (n>=51)  ? arr[51]  : "."
        dbsnp_rs     = (n>=127) ? arr[127] : "."
        drugbank     = (n>=91)  ? arr[91]  : "."
        go_bio       = (n>=94)  ? arr[94]  : "."

        print $1,$2,$3,$4,$5,$6,$7,$8,$9,\
              gene,var_class,sec_class,var_type,\
              cdna_chg,protein_chg,\
              clinvar_sig,clinvar_dis,clinvar_rs,\
              cosmic_mut,cosmic_tot,cosmic_tiss,\
              dbsnp_rs,drugbank,go_bio
    }' >> "$TSV_OUT"

    if [ $? -ne 0 ]; then
        echo "ERROR: TSV extraction failed for P-${i}"
        return 1
    fi

    # Summary
    TOTAL=$(tail -n +2 "$TSV_OUT" | wc -l)
    WITH_RS=$(tail -n +2 "$TSV_OUT" | awk -F'\t' '$22!="." && $22!="" {c++} END{print c+0}')
    WITHOUT_RS=$(tail -n +2 "$TSV_OUT" | awk -F'\t' '$22=="." || $22=="" {c++} END{print c+0}')

    echo ">>> Sample P-${i} final summary:"
    echo "    Total PASS variants:  $TOTAL"
    echo "    Known (rs ID):        $WITH_RS"
    echo "    Novel:                $WITHOUT_RS"


    echo ">>> Completed Sample P-${i} at $(date)"
    echo ">>> Output: $TSV_OUT"
    echo ""
}

export -f run_bestpractices_pipeline
export -f get_bam_paths
