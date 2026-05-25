#!/usr/bin/env bash
#SBATCH --job-name=Admixture
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/resources/RNA/Admixture_job_%j.out
#SBATCH --mem=48G
#SBATCH -c 12

#srun -c 12 --mem=48G --pty -J Admixture bash

cd /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

# Load settings
source ./00_Settings.sh

# Load necessary modules
module load bcftools
module load htslib
module load plink
module load R

temp_dir="/lustre/home/wallbp/tmp"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# Admixture_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_RNA_GATK_0_0"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
# Admixture_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_WGS_GATK_0_0"

HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
Admixture_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_RNA_GATK_resources_test"

# HaplotypeCaller_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
# Admixture_data_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/Admixture_WGS_GATK_resources_test"



#in_1k="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/merged_filtered_1000G.vcf.gz"
#in_1k="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/merged_filtered_1000G.minAFgt5p.vcf.gz"
in_1k="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05.filtered_common_snps.vcf.gz"

mkdir -p "$Admixture_data_dir"

bash /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/admixture_pipeline.sh \
  -t "$in_1k" \
  -i "$HaplotypeCaller_data_dir" \
  -g "/lustre/home/juicer/ref_genomes/jax_refs/human/GRCh38/genome/sequence/gencode/v39/GRCh38.primary_assembly.genome.fa" \
  -a "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/admixture" \
  -r "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/PLOTTING_admixture_results.R" \
  -m "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/20130606_g1k_3202_samples_ped_population.txt" \
  -o "$Admixture_data_dir" \
  -s 42 \
  -p 2 \
  -d "$temp_dir" \
  -n
