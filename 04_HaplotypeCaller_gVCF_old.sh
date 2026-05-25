#!/usr/bin/env bash

# Load settings
source ./00_Settings.sh

for bam in `find ${BaseRecalibrator_data_dir} -type f -name "*.bam"`; do
	
	neutral_base=$(basename "$bam" .bqsr.split.marked.sorted.bam)

    if [ ! -f "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.sorted.g.vcf.gz" ]; then

        echo "$(date '+%Y-%m-%d %H:%M:%S'): Submitting job for $bam"
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=HaplotypeCaller_gVCF_${neutral_base}
#SBATCH --output=./logs/HaplotypeCaller_gVCF_${neutral_base}_%j.out
#SBATCH --error=./logs/HaplotypeCaller_gVCF_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --mem=64G

echo "$(date '+%Y-%m-%d %H:%M:%S'): HaplotypeCaller (gVCF): $bam"
if [ ! -f "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz" ]; then
    gatk HaplotypeCaller \
        --java-options "-XX:ParallelGCThreads=8 -Xmx64G" \
        --input "$bam" \
        --output "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz" \
        --reference $HUMAN_GENOME_REF \
        --adaptive-pruning true \
        --standard-min-confidence-threshold-for-calling 20.0 \
        --emit-ref-confidence GVCF
fi

bcftools sort \
    "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz" \
    --max-mem "64G" \
    --output "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.sorted.g.vcf.gz" \
    --output-type z \
    --write-index=tbi

rm "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz"
rm "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz.tbi"

EOF
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Skipping already processed file: $bam"
        rm "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz"
        rm "$HaplotypeCaller_gvcf_data_dir/${neutral_base}.unsorted.g.vcf.gz.tbi"
    fi
done
