#!bin/bash

module load bcftools

# in_dir_1="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
# in_dir_2="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/compressed"
# out_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_RNA.csv"
# regex='s/_.*|\.vcf.*|\.mpileup.*//g'

# in_dir_1="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
# in_dir_2="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/compressed"
# out_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_WGS.csv"
# regex='s/_.*|\.vcf.*|\.mpileup.*//g'

# in_dir_1="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/HaplotypeCaller"
# in_dir_2="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_RNA_GATK/results/compressed_unfiltered"
# out_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_RNA.unfiltered.csv"
# regex='s/_.*|\.vcf.*|\.mpileup.*//g'

in_dir_1="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF"
in_dir_2="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/compressed_unfiltered"
out_csv="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/GATK_vs_JAX_WGS.unfiltered.csv"
regex='s/_.*|\.vcf.*|\.mpileup.*//g'

declare -A vcf_files_1
for file in $(find "$in_dir_1" -name "*.vcf.gz"); do
  filename=$(basename "$file")

  # RegEx for naming
  key=$(echo "$filename" | sed -E "$regex")
  
  # assign to array
  vcf_files_1["$key"]="$file"
done

declare -A vcf_files_2
for file in $(find "$in_dir_2" -name "*.vcf.gz"); do
  filename=$(basename "$file")

  # RegEx for naming
  key=$(echo "$filename" | sed -E "$regex")
  
  # assign to array
  vcf_files_2["$key"]="$file"
done

# Validate
# for key in "${!vcf_files_1[@]}"; do   echo "$key -> ${vcf_files_1[$key]}"; done
# for key in "${!vcf_files_2[@]}"; do   echo "$key -> ${vcf_files_2[$key]}"; done

# Write header
echo "Sample,Type,Total Records,Total Variants,SNPs,INDELs,MNPs,Other" > $out_csv

declare -A isec_names
for sample in ${!vcf_files_1[@]}; do

    vcf_path_1=${vcf_files_1[$sample]}
    
    # Exit if sample not in in_dir_2
    if [[ ! -v vcf_files_1["$sample"] ]]; then
        echo "Warning - Cannot find ${key} in ${in_dir_2}... skipping..."
        continue
    fi

    vcf_path_2="${vcf_files_2[$sample]}"

    tmp_dir=$(mktemp -d)

    bcftools isec -p "$tmp_dir" "$vcf_path_1" "$vcf_path_2"

    # Count variants
    # 0000.vcf unique_1
    # 0001.vcf unique_2
    # 0002.vcf shared_1
    # 0003.vcf shared_2
    isec_names=(
        [unique_1]="${tmp_dir}/0000.vcf"
        [unique_2]="${tmp_dir}/0001.vcf"
        [shared_1]="${tmp_dir}/0002.vcf"
        [shared_2]="${tmp_dir}/0003.vcf"
    )

    for key in ${!isec_names[@]}; do

        echo "$(date +"%Y-%m-%d %H:%M:%S") Reading ${key} for ${sample}..."

        vcf=${isec_names[$key]}
        TOTAL_RECORDS=$(bcftools view "$vcf" -H | wc -l)
        TOTAL_VARIANTS=$(bcftools view -v snps,indels,mnps,other "$vcf" -H | wc -l)
        SNPS=$(bcftools view -v snps "$vcf" -H | wc -l)
        INDELS=$(bcftools view -v indels "$vcf" -H | wc -l)
        MNPS=$(bcftools view -v mnps "$vcf" -H | wc -l)
        OTHER=$(bcftools view -v other,ref,bnd "$vcf" -H | wc -l)
        
        echo "$sample,$key,$TOTAL_RECORDS,$TOTAL_VARIANTS,$SNPS,$INDELS,$MNPS,$OTHER" >> $out_csv

    done

    rm -r "$tmp_dir"
done

#vcf_path_1=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/JAX_WGS/results/compressed/VCU-LC-099.mpileup.called.filtered.annotated.vcf.vcf.gz
#vcf_path_2=/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered/VCU-LC-099_1795_CKDO240002088-1A.vcf.gz
