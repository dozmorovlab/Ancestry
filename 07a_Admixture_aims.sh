#!bin/bash

module load bcftools

in_vcf="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/Admixture.1k.common_snps.vcf.gz"
out_bed="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models/Admixture.bed"

bcftools norm -m -any "$in_vcf" -Ob | \
    bcftools query -f '%CHROM\t%POS0\t%REF\n' - | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $2 + length($3)}' > "$out_bed"
