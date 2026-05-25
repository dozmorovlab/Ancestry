#!/usr/bin/env bash

# Load settings
cd "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev"
source ./00_Settings.sh

for bam in `find ${BaseRecalibrator_data_dir} -type f -name "*.bam"`; do

    neutral_base=$(basename "$bam" .bqsr.split.marked.sorted.bam)

    if [ ! -f "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.vcf.gz.tbi" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Submitting job for $bam"
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=HaplotypeCaller_gVCF_${neutral_base}
#SBATCH --output=./logs/HaplotypeCaller_gVCF_${neutral_base}_%j.out
#SBATCH --error=./logs/HaplotypeCaller_gVCF_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=32G

source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

gatk HaplotypeCaller \
    --java-options "-XX:ParallelGCThreads=4 -Xmx32G" \
    --input "$bam" \
    --output "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.vcf.gz" \
    --reference $HUMAN_GENOME_REF \
    --adaptive-pruning true \
    --standard-min-confidence-threshold-for-calling 20.0 \
    -ERC GVCF

EOF
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Skipping already processed file: $bam"
    fi
done
