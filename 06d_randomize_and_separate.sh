#!/bin/bash
#SBATCH --job-name=rand_sep
#SBATCH --output=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/rand_sep.out
#SBATCH --error=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/rand_sep.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module load bcftools

MAX_JOBS=100

wait_for_available_slots() {
    while [ "$(squeue -u $USER | grep -cE ' R | PD ')" -ge "$MAX_JOBS" ]; do
        sleep 10
    done
}

in_overlap="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_1k_RNA_filtered.vcf.gz"

out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_RNA_10_select52"
mkdir -p $out_dir

logs="${out_dir}/logs"
mkdir -p $logs

declare -A out_n=(
    ["10"]=10
    ["100"]=100
    ["1k"]=1000
    ["10k"]=10000
    ["100k"]=100000
)

select_1k=(
    "HG00264" "HG00237" "HG00362" "HG00336" "HG00674" "HG00410" "HG01070"
    "HG01167" "HG02186" "HG01816" "HG01383" "HG01435" "HG01705" "HG01670"
    "HG02265" "HG02292" "HG02660" "HG02652" "HG01848" "HG01872" "HG02332"
    "HG01883" "HG02634" "HG03049" "HG03372" "HG03367" "HG04140" "HG03937"
    "HG03091" "HG03376" "HG03684" "HG03851" "HG03742" "HG03977" "NA12045"
    "NA12156" "NA19146" "NA19206" "NA18542" "NA18749" "NA18981" "NA18953"
    "NA19378" "NA19374" "NA20314" "NA19681" "NA19670" "NA20281" "NA20542"
    "NA20509" "NA21113" "NA21110"
)

#rand_seed=42
base="$(basename "$in_overlap" .vcf.gz)"

    for n_name in "${!out_n[@]}"; do
        n=${out_n[$n_name]}

        positions_file="${out_dir}/${n_name}_snps/random_positions.seed-${rand_seed}.n-${n}.txt"
        mkdir -p $(dirname "$positions_file") 

        if [ ! -e $positions_file ]; then
        
            bcftools view \
                -v snps \
                -H \
                "${in_overlap}" | \
            LC_ALL=C shuf \
                --random-source=<(yes "${rand_seed}") \
                -n "${n}" | \
            cut \
                -f1,2 | \
            awk \
                '{print $1 "\t" $2}' \
            > "$positions_file"

        fi

        for sample_1k_name in $(bcftools query -l "${in_overlap}"); do

            # Look in select samples, skip otherwise
            found=0
            for item in "${select_1k[@]}"; do
                if [[ "$item" == "$sample_1k_name" ]]; then
                    found=1
                    break
                fi
            done
            if [[ $found -eq 0 ]]; then
                continue
            fi

            out_file_sample="${out_dir}/${n_name}_snps/${rand_seed}/${sample_1k_name}.${base}.seed-${rand_seed}.n-${n}.vcf.gz"
            mkdir -p "$(dirname "$out_file_sample")"

            out_base="$(basename "$out_file_sample")"
            sample_script="${logs}/${n_name}_snps/${rand_seed}_seed/${out_base}.slurm"
            mkdir -p "$(dirname "$sample_script")"

            if [ ! -e $out_file_sample ]; then
            
                cat > "$sample_script" <<EOF
#!/bin/bash
#SBATCH --job-name=${out_base}
#SBATCH --output=${logs}/${n_name}_snps/${rand_seed}_seed/${out_base}.out
#SBATCH --error=${logs}/${n_name}_snps/${rand_seed}_seed/${out_base}.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

module load bcftools

bcftools view \
    -R "${positions_file}" \
    -v snps \
    -c 1 \
    -s "${sample_1k_name}" \
    -Oz \
    --write-index=tbi \
    -o "${out_file_sample}" \
    "${in_overlap}"
EOF
                wait_for_available_slots
                sbatch "$sample_script"

            fi
        done
    done
done
