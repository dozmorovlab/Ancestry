#!bin/bash

module load bcftools

# Directory containing VCF files
VCF_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_RNA_10_select52"
OUTFILE="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Ancestry_Results/Data/overlap_chrs/random_RNA_10_select52.counts.tsv"

# Write header
echo -e "dir1\tdir2\tfile\tsnp_count" > "$OUTFILE"

# Find all .vcf and .vcf.gz files recursively
find "$VCF_DIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | while read -r vcf; do
    # Get directory path
    dirpath=$(dirname "$vcf")

    # Split into array
    IFS='/' read -r -a parts <<< "$dirpath"

    # Get last two directories (safe with NA if not enough levels)
    n=${#parts[@]}
    dir1="${parts[$((n-2))]:-NA}"
    dir2="${parts[$((n-1))]:-NA}"
    file=$(basename "$vcf")

    # Count SNPs
    count=$(bcftools view -v snps -H "$vcf" | wc -l)
    #count=0

    # Append results
    echo -e "$dir1\t$dir2\t$file\t$count" >> "$OUTFILE"
done

echo "Results written to $OUTFILE"
