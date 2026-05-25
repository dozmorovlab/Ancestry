#!/usr/bin/env bash
#SBATCH --job-name=BQSR_03
#SBATCH --output=./03_BQSR_output.txt
#SBATCH --error=./03_BQSR_errors.txt
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G

# Load settings
source ./00_Settings.sh

#eval known_sites="$known_sites"
known_sites_str=""
for vcf in "${known_sites[@]}"; do
    known_sites_str+="--known-sites ${vcf} "
done

for vcf in "${known_sites[@]}"; do
	# Index known sites once
	if [ ! -f "$vcf.tbi" ]; then
		echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing: $vcf"
		gatk IndexFeatureFile \
			--input $vcf
	fi
done

SplitNCigarReads_bam_dir="${data_dir}/RG"
for bam in `find ${SplitNCigarReads_bam_dir} -type f -name "*.bam"`; do

	echo "$(date '+%Y-%m-%d %H:%M:%S'): Base Quality Recalibration: $bam"

	neutral_base=$(basename "$bam" .split.marked.sorted.bam)

	sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=BQSR_${neutral_base}
#SBATCH --output=./logs/BQSR_${neutral_base}_%j.out
#SBATCH --error=./logs/BQSR_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=4G

echo "\$(date '+%Y-%m-%d %H:%M:%S'): Splitting N cigar reads: $bam"

# BaseRecalibrator
# https://gatk.broadinstitute.org/hc/en-us/articles/13832708374939-BaseRecalibrator
# These Read Filters are automatically applied to the data by the Engine before processing by BaseRecalibrator.
#	NotSecondaryAlignmentReadFilter
#	PassesVendorQualityCheckReadFilter
#	MappedReadFilter
#	MappingQualityAvailableReadFilter
#	NotDuplicateReadFilter
#	MappingQualityNotZeroReadFilter
#	WellformedReadFilter

echo "$(date '+%Y-%m-%d %H:%M:%S'): BaseRecalibrator: $bam"

gatk BaseRecalibrator \
	--input $bam \
	--output "$BaseRecalibrator_data_dir/$neutral_base.BaseRecalibration.table" \
	--reference $HUMAN_GENOME_REF \
	$known_sites_str

echo "$(date '+%Y-%m-%d %H:%M:%S'): AnalyzeCovariates: $bam"

# AnalyzeCovariates
# https://gatk.broadinstitute.org/hc/en-us/articles/13832627271963-AnalyzeCovariates
gatk AnalyzeCovariates \
	-bqsr "$BaseRecalibrator_data_dir/$neutral_base.BaseRecalibration.table" \
	-plots "$BaseRecalibrator_report/$neutral_base.AnalyzeCovariates.pdf"

# ApplyBQSR
# https://gatk.broadinstitute.org/hc/en-us/articles/13832692459163-ApplyBQSR
# This Read Filter is automatically applied to the data by the Engine before processing by ApplyBQSR.
#	WellformedReadFilter

echo "$(date '+%Y-%m-%d %H:%M:%S'): ApplyBQSR: $bam"

gatk ApplyBQSR \
	--input $bam \
	--output "$BaseRecalibrator_data_dir/$neutral_base.bqsr.split.marked.sorted.bam" \
	--bqsr-recal-file "$BaseRecalibrator_data_dir/$neutral_base.BaseRecalibration.table" \
	--reference $HUMAN_GENOME_REF

EOF
done