#!/bin/bash

# Directories
#in_dir="/lustre/home/harrell_lab/WGS/sample_VCFs_c1"
#out_dir="/lustre/home/harrell_lab/WGS/filtered_VCFs"
#log_dir="/lustre/home/harrell_lab/WGS/logs"

#in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/HaplotypeCaller"
#out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/HaplotypeCaller_filtered"
#log_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_filter_logs"

#in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF"
#out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#log_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered_logs"

in_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/HaplotypeCaller"
out_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
log_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered_logs"

mkdir -p "$log_dir" "$out_dir"  # Create directories for logs and output files

for vcf in "$in_dir"/*.vcf.gz; do
    filename=$(basename "$vcf" .vcf.gz)
    echo "Submitting ${filename}"

    # Create a temporary file for the job script
    script=$(mktemp)

    # Write the job script template without variable expansion.
    # Using single quotes around EOF prevents expansion.
    cat <<'EOF' > "$script"
#!/bin/bash
#SBATCH --job-name=FILTER_JOB_NAME
#SBATCH --output=FILTER_LOGDIR/FILENAME.out
#SBATCH --error=FILTER_LOGDIR/FILENAME.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load bcftools

# Function to count SNPs and append data to a TSV file using bcftools
count_snps() {
    local vcf_file="$1"
    local tsv_file="$2"
    local comment="$3"
    local last_count="$4"
    local current_snp_count
    current_snp_count=$(bcftools view "$vcf_file" -H | wc -l)
    local lost_snps=$((current_snp_count - last_count))
    echo -e "$current_snp_count\t$lost_snps\t$comment" >> "$tsv_file"
    echo "$current_snp_count"
}

# Define variables for file paths (placeholders to be replaced)
OUTPUT_VCF="OUTPUT_VCF_PLACEHOLDER"
TSV_FILE_OUT="TSV_FILE_OUT_PLACEHOLDER"
TEMP_VCF_IN="TEMP_VCF_IN_PLACEHOLDER"
TEMP_VCF_OUT="TEMP_VCF_OUT_PLACEHOLDER"
INPUT_VCF="INPUT_VCF_PLACEHOLDER"

# Start processing the VCF
bcftools view "$INPUT_VCF" > "$TEMP_VCF_IN"

last_count=0 

# Define individual arrays for each filter command.
declare -a filter1=( view )
declare -a filter2=( view -m2 -M2 -v snps )
#declare -a filter3=( view -i 'FILTER="PASS"' )
declare -a filter4=( view -e 'GT="./."' )
declare -a filter5=( filter -i 'QUAL>30' )
declare -a filter6=( filter -i 'INFO/DP>10' )
declare -a filter7=( filter -i 'INFO/MQ>30' )
declare -a filter8=( filter -i 'FORMAT/DP>10' )
declare -a filter9=( filter -i 'FORMAT/GQ>20' )
declare -a filter10=( filter -e 'FS>30' )
declare -a filter11=( filter -e 'SOR>3.0' )
declare -a filter12=( filter -e 'GT="0/1" & (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1])) < 0.3' )
declare -a filter13=( filter -e 'GT="0/1" & (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1])) > 0.7' )
declare -a filter14=( filter -e 'GT="1/1" & FORMAT/AD[0:0] / FORMAT/DP > 0.05' )
declare -a filter15=( filter -e 'GT="0/0" & FORMAT/AD[0:1] / FORMAT/DP > 0.05' )

# Create an array containing the names of each command array.
#filter_cmds=( filter1 filter2 filter3 filter4 filter5 filter6 filter7 filter8 filter9 filter10 filter11 filter12 filter13 filter14 filter15 )
filter_cmds=( filter1 filter2 filter4 filter5 filter6 filter7 filter8 filter9 filter10 filter11 filter12 filter13 filter14 filter15 )


for cmd_name in "${filter_cmds[@]}"; do
    # Use eval and indirect expansion to load the command array into "cmd".
    eval "cmd=( \"\${${cmd_name}[@]}\" )"
    bcftools "${cmd[@]}" "$TEMP_VCF_IN" > "$TEMP_VCF_OUT"
    last_count=$(count_snps "$TEMP_VCF_OUT" "$TSV_FILE_OUT" "bcftools ${cmd[*]}" "$last_count")
    mv "$TEMP_VCF_OUT" "$TEMP_VCF_IN"
done

bcftools view -Oz -o "$OUTPUT_VCF" "$TEMP_VCF_IN"
rm "$TEMP_VCF_IN" "$TEMP_VCF_OUT"

echo "Processing completed. Output VCF saved to $OUTPUT_VCF"
EOF

    # Now replace the placeholders with actual values.
    sed -i "s|FILTER_JOB_NAME|filter_${filename}|g" "$script"
    sed -i "s|FILTER_LOGDIR|${log_dir}|g" "$script"
    sed -i "s|FILENAME|${filename}|g" "$script"
    sed -i "s|OUTPUT_VCF_PLACEHOLDER|${out_dir}/${filename}.vcf.gz|g" "$script"
    sed -i "s|TSV_FILE_OUT_PLACEHOLDER|${out_dir}/${filename}_snps.tsv|g" "$script"
    sed -i "s|TEMP_VCF_IN_PLACEHOLDER|${out_dir}/${filename}_temp_in.vcf|g" "$script"
    sed -i "s|TEMP_VCF_OUT_PLACEHOLDER|${out_dir}/${filename}_temp_out.vcf|g" "$script"
    sed -i "s|INPUT_VCF_PLACEHOLDER|${vcf}|g" "$script"

    # Submit the job script
    sbatch "$script"

    # Optionally remove the temporary job script file
    rm "$script"
done
