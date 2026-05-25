#!/bin/bash

home_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev"

separate_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_RNA_10_select52"

out_pred_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/predicted_RNA_10_select52"
mkdir -p $out_pred_dir

logs="${out_pred_dir}/logs"
mkdir -p $logs

declare -A out_n=(
    ["10"]=10
    ["100"]=100
    ["1k"]=1000
    ["10k"]=10000
    ["100k"]=100000
)

for rand_seed in {1..10}; do
    for n_name in "${!out_n[@]}"; do
        n=${out_n[$n_name]}

        dir="${separate_dir}/${n_name}_snps/${rand_seed}"
        dir_name="$n_name"

        [ -d "$dir" ] || continue

        script_dir="${logs}/job_scripts/${dir_name}_snps/${rand_seed}"
            mkdir -p "$script_dir"

        # RAIDS
        raids_job_script="${script_dir}/predict_raids_${dir_name}_snps_${rand_seed}.sh"

        raids_cpus=2
        raids_mem=32G

        cat <<EOL > "$raids_job_script"
#!/bin/bash
#SBATCH --job-name=predict_raids_${dir_name}_snps_${rand_seed}
#SBATCH --output=${logs}/predict_raids_${dir_name}_snps_${rand_seed}.out
#SBATCH --error=${logs}/predict_raids_${dir_name}_snps_${rand_seed}.err
#SBATCH -c ${raids_cpus}
#SBATCH --mem=${raids_mem}

# Load required modules (example for bcftools, adjust as needed)
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

cd "$home_dir"
source "00_Settings.sh"

HaplotypeCaller_data_dir="$separate_dir/${dir_name}_snps/${rand_seed}"

RAIDS_data_dir="$out_pred_dir/${dir_name}_snps/${rand_seed}/RAIDS"
mkdir -p \$RAIDS_data_dir

Rscript "$home_dir/run_RAIDS.R" \
    --input_pop_gds "\$matGeno1000g" \
    --input_annot_gds "\$matAnnot1000g" \
    --input_dir "\$HaplotypeCaller_data_dir" \
    --output_gds "\$RAIDS_data_dir" \
    --output_results "\$RAIDS_data_dir" \
    --threads "$raids_cpus" \
    --n_profiles "\$number_of_profiles"

EOL

        # EthSEQ
        ethseq_job_script="${script_dir}/predict_ethseq_${dir_name}_snps_${rand_seed}.sh"

        ethseq_cpus=2
        ethseq_mem=32G

        cat <<EOL > "$ethseq_job_script"
#!/bin/bash
#SBATCH --job-name=predict_ethseq_${dir_name}_snps_${rand_seed}
#SBATCH --output=${logs}/predict_ethseq_${dir_name}_snps_${rand_seed}.out
#SBATCH --error=${logs}/predict_ethseq_${dir_name}_snps_${rand_seed}.err
#SBATCH -c ${ethseq_cpus}
#SBATCH --mem=${ethseq_mem}

# Load required modules (example for bcftools, adjust as needed)
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

cd "$home_dir"
source "00_Settings.sh"

HaplotypeCaller_data_dir="$separate_dir/${dir_name}_snps/${rand_seed}"

EthSEQ_data_dir="$out_pred_dir/${dir_name}_snps/${rand_seed}/EthSEQ"
mkdir -p "\$EthSEQ_data_dir"

Rscript "$home_dir/run_EthSEQ.R" \
    --input_dir "\$HaplotypeCaller_data_dir" \
    --output_results "\$EthSEQ_data_dir" \
    --cores "$ethseq_cpus" \
    --pops "All"

EOL

        # gnomAD
        gnomad_job_script="${script_dir}/predict_gnomad_${dir_name}_snps_${rand_seed}.sh"

        gnomad_cpus=2
        gnomad_mem=64G

        cat <<EOL > "$gnomad_job_script"
#!/bin/bash
#SBATCH --job-name=predict_gnomad_${dir_name}_snps_${rand_seed}
#SBATCH --output=${logs}/predict_gnomad_${dir_name}_snps_${rand_seed}.out
#SBATCH --error=${logs}/predict_gnomad_${dir_name}_snps_${rand_seed}.err
#SBATCH -c ${gnomad_cpus}
#SBATCH --mem=${gnomad_mem}

cd "$home_dir"
source "00_Settings.sh"

# Load required modules (example for bcftools, adjust as needed)
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/gnomad

HaplotypeCaller_data_dir="$separate_dir/${dir_name}_snps/${rand_seed}"
gnomAD_data_dir="$out_pred_dir/${dir_name}_snps/${rand_seed}/gnomAD"
mkdir -p \$gnomAD_data_dir

python3 run_gnomAD.py \
    --vcf_dir "\$HaplotypeCaller_data_dir" \
    --sequence_report "\$HUMAN_SEQ_REPORT" \
    --loadings_path "\$loadings_path" \
    --rf_model_path "\$rf_model_path" \
    --hgdp_1kg_mt "\$hgdp_1kg_mt_path" \
    --out_hgdp_1kg_pred_path "\${gnomAD_data_dir}/gnomAD_hgdp_1kg_pred.rna_filtered.csv" \
    --out_sample_pred_path "\${gnomAD_data_dir}/gnomAD_sample_pred.rna_filtered.csv" \
    --temp_dir "\${gnomAD_data_dir}/temp" \
    --num_pcs 20 \
    --min_prob 0.0 \
    --file_suffix ".vcf.gz" \
    --mem "$gnomad_mem" \
    --threads "$gnomad_cpus"

EOL

        # Admixture
        admixture_job_script="${script_dir}/predict_admixture_${dir_name}_snps_${rand_seed}.sh"

        admixture_cpus=12
        admixture_mem=24G

        cat <<EOL > "$admixture_job_script"
#!/bin/bash
#SBATCH --job-name=predict_admixture_${dir_name}_snps_${rand_seed}
#SBATCH --output=${logs}/predict_admixture_${dir_name}_snps_${rand_seed}.out
#SBATCH --error=${logs}/predict_admixture_${dir_name}_snps_${rand_seed}.err
#SBATCH -c ${admixture_cpus}
#SBATCH --mem=${admixture_mem}

cd "$home_dir"
source "00_Settings.sh"

# Load necessary modules
module load bcftools
module load htslib
module load plink
module load R

HaplotypeCaller_data_dir="${separate_dir}/${dir_name}_snps/${rand_seed}"
Admixture_data_dir="${out_pred_dir}/${dir_name}_snps/${rand_seed}/Admixture"

mkdir -p "\$Admixture_data_dir"

bash /lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/admixture_pipeline.sh \
-t "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/merged_filtered_1000G.vcf.gz" \
-i "\$HaplotypeCaller_data_dir" \
-g "/lustre/home/juicer/ref_genomes/jax_refs/human/GRCh38/genome/sequence/gencode/v39/GRCh38.primary_assembly.genome.fa" \
-a "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/admixture" \
-r "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/PLOTTING_admixture_results.R" \
-m "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture/20130606_g1k_3202_samples_ped_population.txt" \
-o "\$Admixture_data_dir" \
-s 42 \
-p "${admixture_cpus}"

EOL

        # Aeon
        aeon_job_script="${script_dir}/predict_aeon_${dir_name}_snps_${rand_seed}.sh"

        aeon_cpus=2
        aeon_mem=64G

        cat <<EOL > "$aeon_job_script"
#!/bin/bash
#SBATCH --job-name=predict_aeon_${dir_name}_snps_${rand_seed}
#SBATCH --output=${logs}/predict_aeon_${dir_name}_snps_${rand_seed}.out
#SBATCH --error=${logs}/predict_aeon_${dir_name}_snps_${rand_seed}.err
#SBATCH -c ${aeon_cpus}
#SBATCH --mem=${aeon_mem}

cd "$home_dir"
source "00_Settings.sh"

source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate aeon
module load bcftools

HaplotypeCaller_data_dir="${separate_dir}/${dir_name}_snps/${rand_seed}"
Aeon_data_dir="${out_pred_dir}/${dir_name}_snps/${rand_seed}/Aeon"

mkdir -p \$Aeon_data_dir

echo "\$(date '+%Y-%m-%d %H:%M:%S'): Merging VCFs..."
TEMP_DIR=\$(mktemp -d)
bcftools merge "\$HaplotypeCaller_data_dir"/*.vcf.gz -W -Oz -o "\${TEMP_DIR}/merged.vcf.gz"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Running Aeon..."
cd "\$Aeon_data_dir"
aeon \
    "\${TEMP_DIR}/merged.vcf.gz" \
    -o "\$Aeon_data_dir/\$(basename \${Aeon_data_dir})"

EOL

        # JAX - TBD

        # Submit the jobs using sbatch
        sbatch "$raids_job_script"
        #sbatch "$ethseq_job_script"
        #sbatch "$gnomad_job_script"
        #sbatch "$admixture_job_script"
        #sbatch "$aeon_job_script"

    done
done

