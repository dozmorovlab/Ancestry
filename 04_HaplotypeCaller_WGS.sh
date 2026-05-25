#!/usr/bin/env bash

# Load settings
source ./00_Settings.sh
BaseRecalibrator_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS"
HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF"

mkdir -p $HaplotypeCaller_data_dir

for bam in "${BaseRecalibrator_data_dir}/"*.bam; do
    echo $bam
    neutral_base=$(basename "$bam" .merged.sorted.bqsr.bam)

    if [ ! -f "$HaplotypeCaller_data_dir/${neutral_base}.vcf.gz" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Submitting job for $bam"
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=HaplotypeCaller_${neutral_base}
#SBATCH --output=./logs/HaplotypeCaller_${neutral_base}_%j.out
#SBATCH --error=./logs/HaplotypeCaller_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=16G

gatk HaplotypeCaller \
    --java-options "-XX:ParallelGCThreads=4 -Xmx16G" \
    --input "$bam" \
    --output "$HaplotypeCaller_data_dir/${neutral_base}.vcf.gz" \
    --reference $HUMAN_GENOME_REF \
    --adaptive-pruning true \
    --standard-min-confidence-threshold-for-calling 20.0

EOF
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Skipping already processed file: $bam"
    fi
done