#!/bin/bash

# Directories
in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
#in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"

module load bcftools

for vcf in "$in_dir"/*.vcf.gz; do
  
  echo "Indexing $(basename $vcf)..."
  bcftools index $vcf

done
