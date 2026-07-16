#All commands for running freeC and getting the genes in a single file should be ran one by one
#to run freec
mamba activate freeC_env
./New	

#after freeC to get genes

#unzip annotation
gunzip gencode.v21.chr_patch_hapl_scaff.annotation.gtf.gz 

#get data as a bed file for hg38
awk '$3=="gene" {
    match($0, /gene_name "([^"]+)"/, arr);
    print $1"\t"$4-1"\t"$5"\t"arr[1]
}' gencode.v21.chr_patch_hapl_scaff.annotation.gtf > genes.bed

head genes.bed

#if you use the alternate annotation use this - I used this 
gunzip gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz 
awk '$3=="gene" {                                                              
	match($0, /gene_name "([^"]+)"/, arr);
	print $1"\t"$4-1"\t"$5"\t"arr[1]
}' gencode.v49.chr_patch_hapl_scaff.annotation.gtf > genes_49.bed
 head genes_49.bed

#save the .bam_CNVs as a cnv.bed
#to do this in a loop for all files in the parent folder (starting with T) and ending with .bam_CNVs
# Run from parent folder with T* directories
for dir in T*; do
    cd "$dir"
    
# Find the .bam_CNVs file and convert to standard BED format
    if ls *.bam_CNVs >/dev/null 2>&1; then
        awk 'BEGIN{OFS="\t"}
             {
                chr=$1; start=$2; end=$3; name=$4; 
                # Add chr prefix only if missing
                if (chr !~ /^chr/) { chr="chr"chr }
                print chr, start, end, name
             }' *.bam_CNVs > cnv.bed
        
        echo "Created cnv.bed in $dir (standard BED4 format)"
    fi
    cd ..
done






# Copy genes_49.bed ONLY to directories that contain cnv.bed files
for dir in T*; do
    if [ -d "$dir" ] && [ -f "$dir/cnv.bed" ]; then
        cp /home/aarti/NII/freec/no_duplicates/genes_49.bed "$dir/"
        echo "Copied genes_49.bed to $dir"
    fi
done


# Annotate CNVs in all T* folders
for dir in T*; do
    if [ -f "$dir/cnv.bed" ] && [ -f "$dir/genes_49.bed" ]; then
        cd "$dir"
        bedtools intersect -a cnv.bed -b genes_49.bed -wa -wb > cnv_genes_49.tsv
        echo "Annotated: $dir"
        cd ..
    fi
done


#to make a single file with all the cnv_genes_49.tsv for all smaples and adding the foldername as a separate column
awk 'BEGIN{OFS="\t"}
FNR==1{
    folder=FILENAME
    sub(/\/cnv_genes_49\.tsv$/, "", folder)
    sub(/^.*\//, "", folder)
}
{print folder, $0}' T*/cnv_genes_49.tsv > combined_cnv_genes_49.tsv


#plotting was done in R


#to view bed graph files in genome browser
#some chr lenghts ar elonger than chr in UCSC..bec of alt chr. So we need to clip it to UCSC size
# Download hg38 chrom.sizes
wget -O hg38.chrom.sizes http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes






#to combine all the ratio.txt files
(
first=1

for f in T*/non_bacterial*.bam_ratio.txt
do
    sample=$(basename $(dirname "$f"))

    if [ $first -eq 1 ]; then
        awk -v s="$sample" '
        BEGIN{OFS="\t"}
        NR==1 {print "sample",$0}
        NR>1  {print s,$0}
        ' "$f"

        first=0
    else
        awk -v s="$sample" '
        BEGIN{OFS="\t"}
        NR>1 {print s,$0}
        ' "$f"
    fi

done
) > merged_ratio.txt




#for scoring
#to make a single file with all the .bam_CNVs for all samples and adding the foldername as a separate column
awk 'BEGIN{OFS="\t"}
FNR==1{
    folder=FILENAME
    sub(/\/[^/]+$/, "", folder)
    sub(/^.*\//, "", folder)
}
{
    print folder, $0
}' T*/non_bacterial*bam_CNVs > combined_non_bacterial_bam_CNVs.tsv


#add headers and chr   if needed
awk 'BEGIN{OFS="\t"}
{
cnv=$1; cnv_start=$2; cnv_end=$3; status=$4; bins=$5; median_ratio=$6; mean_ratio=$7; log2ratio=$8; sample=$9; cnv_id=$10;
# Add chr prefix in column 1

cnv_chr = "chr" cnv
print cnv_chr, cnv_start, cnv_end, status, bins, median_ratio, mean_ratio, log2ratio, sample, cnv_id, Copies

 }' CNV_summary_by_sample.tsv > cnv_filtered_corrected.bed
 
 bedtools intersect -a cnv_filtered_corrected.bed -b genes_49.bed -wa -wb > cnv_scored_annotated.tsv
 
 
 
 #for intgrating cytobands
bedtools intersect -a /PathTO/freec/no_duplicates/results/cnv_gene_collapsed_table.tsv -b cytoBand.txt.gz -wa -wb | cut -f 1-3,8 > cytoband_combined_results.tsv
 

