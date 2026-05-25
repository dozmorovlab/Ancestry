#!/bin/bash

# https://doi.org/10.1158/1055-9965.EPI-18-1132

module load bcftools

threads=2

BED_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models"
#BED_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/papers/bed"

#VCF_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA_VCF_filtered"
#OUT_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/RNA_paper"
#OUT_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/RNA"

VCF_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF_filtered"
#OUT_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/WGS_paper"
OUT_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/results/WGS"

mkdir -p "$OUT_DIR"

for BED in "$BED_DIR"/*.bed; do
    bed_base=$(basename "$BED" .bed)
    echo "$(date) Preparing jobs for ${bed_base}..."

    for vcf in "$VCF_DIR"/*.vcf.gz; do
        vcf_base=$(basename "$vcf" .vcf.gz)
        out="$OUT_DIR/${vcf_base}.${bed_base}.aims.bed"

        if [[ -f "$out" ]]; then
            echo "$(date) Output file $out already exists. Skipping..."
            continue
        fi

        job_script=$(mktemp)
        cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --job-name=${vcf_base}_${bed_base}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem=4G

module load bcftools

if [[ ! -f "$vcf.tbi" ]]; then
    bcftools index -t "$vcf"
fi

bcftools view \
    --threads "$threads" \
    -R "$BED" \
    "$vcf" \
    -Ou | \
bcftools query \
    -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\n' \
    - \
> "$out"
EOF

        sbatch "$job_script"
    done
done
