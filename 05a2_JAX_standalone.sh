#!/bin/bash

srun -c 10 --mem=194G --pty -J JAX bash
cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source ./00_Settings.sh

module load samtools
module load htslib
module load bcftools

source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/nf

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA"
#JAX_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA"
#JAX_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/JAX_samplesheet_RNA.csv"

#HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS"
#JAX_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS"
#JAX_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/JAX_samplesheet_WGS.csv"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/BaseRecalibrator"
JAX_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK"
JAX_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/JAX_samplesheet_RNA_GATK.csv"

mkdir -p "$JAX_data_dir"

work_dir="${JAX_data_dir}/work"
results_dir="${JAX_data_dir}/results"

mkdir -p "$work_dir"
mkdir -p "$results_dir"

# module load nextflow
# Nextflow version 23.10.0 does not match workflow required version: >=24.04.0


#conda activate nf

# nextflow run /lustre/home/harrell_lab/JAX/cs-nf-pipelines \
#     -c /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/apollo_ancestry.config \
#     --workflow ancestry \
#     --csv_input "$JAX_csv" \
#     --pubdir "$JAX_data_dir" \
#     --sample_folder "$HaplotypeCaller_data_dir" \
#     -w $work_dir \
#     --ref_fa /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/References/GRCh38/genome.fa \
#     --genotype_targets /lustre/home/harrell_lab/JAX/ref_data/snp_panel_v2_targets_annotations.snpwt.bed.gz \
#     --snpID_list /lustre/home/harrell_lab/JAX/ref_data/snp_panel_v2.list \
#     --snp_annotations /lustre/home/harrell_lab/JAX/ref_data/snp_panel_v2_targets_annotations.snpwt.bed.gz \
#     --snpweights_panel /lustre/home/harrell_lab/JAX/ref_data/ancestry_panel_v2.snpwt
# For users external to JAX, reference data used by the workflows is available in a Google Cloud bucket, and transfers of data can be made upon request. 
# https://github.com/TheJacksonLaboratory/cs-nf-pipelines/wiki/
# /lustre/home/harrell_lab/JAX/ref_data/get_ref_data.sh

# Removed
# VCU-BC-057,/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/BaseRecalibrator/VCU-BC-057_1000_human.genome.sorted.bqsr.split.marked.sorted.bam