#!/bin/bash

#in_merged_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/RNA_VCF_filtered.merged.vcf.gz"
#out_overlap_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_1k_RNA_filtered.vcf.gz"

in_merged_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/WGS_VCF_filtered.merged.vcf.gz"
out_overlap_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_1k_WGS_filtered.vcf.gz"

k_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/1kgenomes"

#work_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs_RNA"
work_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs_WGS"
mkdir -p "$work_dir"

temp_dir="${work_dir}/temp"
mkdir -p $temp_dir

logs="${work_dir}/logs"
mkdir -p "$logs"

module load bcftools

merged_base=$(basename "$in_merged_file")
for chr_file in $k_dir/*.vcf.gz; do
    chr_base=$(basename $chr_file)
    job_script="${logs}/bcftools_overlap_${chr_base}.slurm"
    
    cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=bcftools_overlap_${chr_base}
#SBATCH --output=${logs}/bcftools_overlap_${chr_base}.out
#SBATCH --error=${logs}/bcftools_overlap_${chr_base}.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=2

set -euo pipefail

bcftools view \\
    --threads 2 \\
    -R "${in_merged_file}" \\
    -O z \\
    -o "${temp_dir}/${merged_base}_${chr_base}.overlap.vcf.gz" \\
    "${chr_file}"
EOF

    # Submit the job and store its ID
    jobid=$(sbatch --parsable "$job_script")
    chr_jobs+=("$jobid")

done

# Create dependency string for merge job
deps=$(IFS=:; echo "${chr_jobs[*]}")
merge_script="${logs}/bcftools_merge.slurm"

cat > "$merge_script" <<EOF
#!/bin/bash
#SBATCH --job-name=bcftools_merge
#SBATCH --output=${logs}/bcftools_merge.out
#SBATCH --error=${logs}/bcftools_merge.err
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

set -euo pipefail

bcftools concat \\
    --regions-overlap 0 \\
    ${temp_dir}/*.overlap.vcf.gz | \\
bcftools sort \\
    --write-index=tbi \\
    -O z \\
    -o "${out_overlap_file}"

EOF

# Submit the merge job after all chr jobs are done
sbatch --dependency=afterok:$deps "$merge_script"

