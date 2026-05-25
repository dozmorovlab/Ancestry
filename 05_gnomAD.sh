#!/usr/bin/env bash
#SBATCH --job-name=gnomAD
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/WGS/gnomAD_job_%j.out
#SBATCH --mem=300G
#SBATCH -c 22

#srun -c 22 --mem=300G --pty -J gnomAD bash
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/gnomad

# Load settings
source ./00_Settings.sh

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA/results/compressed"
#gnomAD_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_RNA"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#gnomAD_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_WGS_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# gnomAD_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_RNA_GATK"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# gnomAD_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_RNA_GATK_resources_test"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
gnomAD_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/gnomAD_WGS_GATK_resources_test"

mkdir -p $gnomAD_data_dir

python3 run_gnomAD.py \
    --vcf_dir "$HaplotypeCaller_data_dir" \
    --sequence_report "$HUMAN_SEQ_REPORT" \
    --loadings_path "$loadings_path" \
    --rf_model_path "$rf_model_path" \
    --hgdp_1kg_mt "$hgdp_1kg_mt_path" \
    --out_hgdp_1kg_pred_path "${gnomAD_data_dir}/gnomAD_hgdp_1kg_pred.csv" \
    --out_sample_pred_path "${gnomAD_data_dir}/gnomAD_sample_pred.csv" \
    --temp_dir "${gnomAD_data_dir}/temp" \
    --num_pcs 20 \
    --min_prob 0.0 \
    --file_suffix ".vcf.gz" \
    --mem 300G \
    --threads 22
