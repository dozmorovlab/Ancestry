#!/bin/bash

home_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev"

separate_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_multi_10"

out_pred_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/predicted_multi_7_10_2025"
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

# List of base directories with BED files
for n_name in "${!out_n[@]}"; do
    n=${out_n[$n_name]}

    dir="${separate_dir}/${n_name}_snps"
    dir_name="$n_name"

    [ -d "$dir" ] || continue

    script_dir="${logs}/job_scripts/${dir_name}_snps"
        mkdir -p "$script_dir"

    # EthSEQ
    ethseq_job_script="${script_dir}/predict_ethseq_${dir_name}_snps.sh"

    ethseq_cpus=2
    ethseq_mem=32G

    cat <<EOL > "$ethseq_job_script"
#!/bin/bash
#SBATCH --job-name=predict_ethseq_${dir_name}_snps
#SBATCH --output=${logs}/predict_ethseq_${dir_name}_snps.out
#SBATCH --error=${logs}/predict_ethseq_${dir_name}_snps.err
#SBATCH -c ${ethseq_cpus}
#SBATCH --mem=${ethseq_mem}

# Load required modules (example for bcftools, adjust as needed)
source /lustre/home/wallbp/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry

cd "$home_dir"
source "00_Settings.sh"

HaplotypeCaller_data_dir="$separate_dir/${dir_name}_snps"

EthSEQ_data_dir="$out_pred_dir/${dir_name}_snps/EthSEQ"
mkdir -p "\$EthSEQ_data_dir"

Rscript "$home_dir/run_EthSEQ.R" \
    --input_dir "\$HaplotypeCaller_data_dir" \
    --output_results "\$EthSEQ_data_dir" \
    --cores "$ethseq_cpus" \
    --pops "All"

EOL

    # Submit the jobs using sbatch
    sbatch "$ethseq_job_script"

done
