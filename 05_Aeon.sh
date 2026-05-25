#!/usr/bin/env bash
#SBATCH --job-name=Aeon
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/WGS/Aeon_job_%j.out
#SBATCH --mem=64G
#SBATCH -c 2

#srun -c 2 --mem=64G --pty -J Aeon bash
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev

# Load settings
#source ./00_Settings.sh

source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate aeon
module load bcftools

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA/results/compressed"
#Aeon_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_RNA"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#Aeon_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_WGS_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# Aeon_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_RNA_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# Aeon_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_RNA_GATK_resources_test"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
Aeon_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Aeon_WGS_GATK_resources_test"

mkdir -p $Aeon_data_dir

echo "$(date '+%Y-%m-%d %H:%M:%S'): Merging VCFs..."
TEMP_DIR=$(mktemp -d)
bcftools merge "$HaplotypeCaller_data_dir"/*.vcf.gz -W -Oz -o "${TEMP_DIR}/merged.vcf.gz"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Running Aeon..."
cd "$Aeon_data_dir"
aeon \
    "${TEMP_DIR}/merged.vcf.gz" \
    -o "$Aeon_data_dir/$(basename ${Aeon_data_dir})"
