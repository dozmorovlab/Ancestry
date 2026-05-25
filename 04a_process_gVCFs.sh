#!/bin/bash
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --output=./logs/GenotypeGVCFs_%j.out
#SBATCH --error=./logs/GenotypeGVCFs_%j.err
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=32G

# Load settings
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source ./00_Settings.sh

source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

echo "$(date '+%Y-%m-%d %H:%M:%S'): Combining gVCFs..."
gatk CombineGVCFs \
    --java-options "-XX:ParallelGCThreads=4 -Xmx32G" \
    --reference "$HUMAN_GENOME_REF" \
    $(printf -- "--variant %s " "${HaplotypeCaller_gvcf_data_dir}/"*) \
    --output "${GenotypeGVCFs_data_dir}/combined.vcf.gz"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Genotyping gVCFs..."
gatk GenotypeGVCFs \
    --java-options "-XX:ParallelGCThreads=4 -Xmx32G" \
    --reference "$HUMAN_GENOME_REF" \
    --variant "${GenotypeGVCFs_data_dir}/combined.vcf.gz" \
    --output "${GenotypeGVCFs_data_dir}/combined.genotyped.vcf.gz"
