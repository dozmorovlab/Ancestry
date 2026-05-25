#!/bin/bash

module load bcftools

# in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/bcftools_annotate"
# out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/compressed"

# in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/bcftools_annotate"
# out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/compressed"

in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/bcftools_call"
out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/compressed_unfiltered"

# in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/bcftools_call"
# out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/compressed_unfiltered"

mkdir -p "$out_dir"

for vcf in "$in_dir"/*.vcf; do
    base=$(basename "$vcf" .vcf)

    echo "$(date) Processing: ${base}..."

    bcftools view \
        "$vcf" \
        -Oz \
        -o "${out_dir}/${base}.vcf.gz" \
        --write-index=tbi
done
