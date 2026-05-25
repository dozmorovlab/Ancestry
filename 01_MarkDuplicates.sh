#!/usr/bin/env bash
#SBATCH --job-name=MarkDuplicates_01
#SBATCH --output=./01_MarkDuplicates.txt
#SBATCH --error=./01_MarkDuplicates_errors.txt
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G

# https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
# https://github.com/nf-core/rnaseq/blob/master/conf/modules.config

# Load settings
source ./00_Settings.sh

max_mem_MB=$(($MAX_MEM_GB * 1024))
current_limit=$(ulimit -n)
max_open_files=$((current_limit - 1000))

for bam in `find "$BAM_DIR" -type f -name "$BAM_REGEX"`; do

    neutral_base=$(basename "$bam" .bam)

    # !!! Assumes input BAM files are sorted, indexed, and properly processed

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Marking duplicates: $bam"
    picard \
        MarkDuplicates \
        --INPUT "${bam}" \
        --OUTPUT "$MarkDuplicates_bam_dir/$neutral_base.marked.sorted.bam" \
        --REFERENCE_SEQUENCE $HUMAN_GENOME_REF \
        --METRICS_FILE "$MarkDuplicates_report/$neutral_base.MarkDuplicates.metrics.txt" \
        -Xmx${MAX_MEM_GB}G \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP "$max_open_files" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE "$OPTICAL_DUPLICATE_PIXEL_DISTANCE" \
        --READ_NAME_REGEX null \
        --REMOVE_DUPLICATES false \
        --COMPRESSION_LEVEL 9
    rm "${neutral_base}_temp.bam"
done