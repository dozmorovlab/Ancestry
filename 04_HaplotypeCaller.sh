#!/usr/bin/env bash

# Load settings
source ./00_Settings.sh

for bam in `find ${BaseRecalibrator_data_dir} -type f -name "*.bam"`; do

    neutral_base=$(basename "$bam" .bqsr.split.marked.sorted.bam)

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
