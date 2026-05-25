#!bin/bash

source ./00_Settings.sh

out_dir=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/RG
log_dir=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/RG/logs

mkdir -p "$out_dir"
mkdir -p "$log_dir"

for bam in $(find ${SplitNCigarReads_bam_dir} -type f -name "*.bam"); do
    neutral_base=$(basename "$bam" .split.marked.sorted.bam)
    
    if [ ! -f "${out_dir}/${neutral_base}.split.marked.sorted.bam" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Submitting job for $bam"
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Add_RG_${neutral_base}
#SBATCH --output=${log_dir}/Add_RG_${neutral_base}_%j.out
#SBATCH --error=${log_dir}/Add_RG_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=4G

echo "\$(date '+%Y-%m-%d %H:%M:%S'): Add_RG: $bam"

picard AddOrReplaceReadGroups \
    I=$bam \
    O=${out_dir}/${neutral_base}.split.marked.sorted.bam \
    RGID=${neutral_base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=${neutral_base}

samtools index \
    ${out_dir}/${neutral_base}.split.marked.sorted.bam
EOF
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Skipping already processed file: $bam"
    fi
done
