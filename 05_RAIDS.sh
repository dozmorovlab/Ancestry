#!/usr/bin/env bash
#SBATCH --job-name=RAIDS
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/RNA/RAIDS_job_%j.out
#SBATCH --mem=71G
#SBATCH -c 6

#srun -c 6 --mem=71G --pty -J RAIDS bash
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

# Load settings
source ./00_Settings.sh

ulimit -n 65535

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA/results/compressed"
#RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_RNA"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_WGS_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RAIDS_leftovers"
# RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RAIDS_leftovers/results"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
#RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_RNA_GATK"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_RNA_GATK_resources_test"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
# RAIDS_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RAIDS_WGS_GATK_resources_test"

mkdir -p $RAIDS_data_dir

MAX_CPUS=6

Rscript run_RAIDS.R \
    --input_pop_gds "$matGeno1000g" \
    --input_annot_gds "$matAnnot1000g" \
    --input_dir "$HaplotypeCaller_data_dir" \
    --output_gds "$RAIDS_data_dir" \
    --output_results "$RAIDS_data_dir" \
    --threads $MAX_CPUS \
    --n_profiles $number_of_profiles
