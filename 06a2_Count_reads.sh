#!/bin/bash

threads=30

# Directory containing BAM files
#bam_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/BaseRecalibrator"
bam_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS"

# Output TSV file
#output_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/RNA_GATK_reads.tsv"
output_file="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results_GATK/Data/WGS_GATK_reads.tsv"

module load samtools

# Write header to the output file
echo -e "file\ttotal_reads" > "$output_file"

# Loop over matching BAM files
for file in "$bam_dir"/*.bam; do

    echo "$(date) Counting reads ${file}..."

    # Count reads using samtools
    count=$(samtools view -@ "$threads" -c "$file")

    # Output the key and sum to TSV
    echo -e "$file\t$count" >> "$output_file"
done
