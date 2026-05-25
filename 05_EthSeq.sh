#!/usr/bin/env bash
#SBATCH --job-name=EthSEQ
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/WGS/EthSEQ_job_%j.out
#SBATCH --mem=120G
#SBATCH -c 8

#srun -c 8 --mem=194G --pty -J EthSEQ bash
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

MAX_CPUS=8

# Load settings
#source ./00_Settings.sh

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA/results/compressed"
#EthSEQ_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_RNA"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#EthSEQ_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_WGS_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# EthSEQ_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_RNA_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# EthSEQ_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_RNA_GATK_resources_test"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
EthSEQ_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/EthSEQ_WGS_GATK_resources_test"

mkdir -p $EthSEQ_data_dir

echo "$(date '+%Y-%m-%d %H:%M:%S'): Running EthSEQ"
Rscript run_EthSEQ.R \
    --input_dir "$HaplotypeCaller_data_dir" \
    --output_results "$EthSEQ_data_dir" \
    --cores "$MAX_CPUS" \
    --pops "All"
