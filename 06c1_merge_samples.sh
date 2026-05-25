#!/bin/bash

#in_vcf_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
#out_merged_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RNA_VCF_filtered.merged.vcf.gz"

in_vcf_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
out_merged_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/WGS_VCF_filtered.merged.vcf.gz"

module load bcftools

# Index all .vcf.gz files if index is missing
for file in "$in_vcf_dir"/*.vcf.gz; do
    if [[ ! -f "$file.tbi" ]]; then
        echo "Indexing $file"
        bcftools index -t "$file"
    fi
done

# Merge all .vcf.gz files
echo "Merging VCFs..."
bcftools merge \
    -Oz \
    --threads 1 \
    --write-index=tbi \
    -o "$out_merged_file" \
    "$in_vcf_dir"/*.vcf.gz

echo "Merged file written to $out_merged_file"
