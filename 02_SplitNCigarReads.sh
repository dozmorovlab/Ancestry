#!/usr/bin/env bash
#SBATCH --job-name=SplitNCigarReads_02            # Job name
#SBATCH --output=./02_SplitNCigarReads.txt        # Standard output log
#SBATCH --error=./02_SplitNCigarReads_errors.txt  # Standard error log
#SBATCH --ntasks=1                                # Number of tasks (processes)
#SBATCH --cpus-per-task=48                        # Number of CPU cores per task
#SBATCH --mem=300G                                # Memory allocation

# Load settings
source ./00_Settings.sh

# Create sequence dictionary once
if [ ! -f "$HUMAN_GENOME_REF.dict" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Creating sequence dictionary: $HUMAN_GENOME_REF"
    gatk CreateSequenceDictionary -R $HUMAN_GENOME_REF
fi

# Index reference once
if [ ! -f "$HUMAN_GENOME_REF.fai" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing reference: $HUMAN_GENOME_REF"
    samtools faidx "$HUMAN_GENOME_REF"
fi

# Create a directory for logs if it doesn't exist
mkdir -p logs

# Iterate over BAM files and check for existing output
for bam in $(find "$MarkDuplicates_bam_dir" -type f -name "*.bam"); do
    neutral_base=$(basename "$bam" .marked.sorted.bam)
    
    # Check if output file already exists
    if [ ! -f "$SplitNCigarReads_bam_dir/$neutral_base.split.marked.sorted.bam" ]; then
        # Submit the SLURM job
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=SplitNCigar_${neutral_base}
#SBATCH --output=./logs/SplitNCigarReads_${neutral_base}_%j.out
#SBATCH --error=./logs/SplitNCigarReads_${neutral_base}_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=4G

echo "\$(date '+%Y-%m-%d %H:%M:%S'): Splitting N cigar reads: $bam"

gatk SplitNCigarReads \
    --input "$bam" \
    --output "$SplitNCigarReads_bam_dir/$neutral_base.split.marked.sorted.bam" \
    --reference "$HUMAN_GENOME_REF"
EOF
    else
        echo "\$(date '+%Y-%m-%d %H:%M:%S'): Skipping already processed file: $bam"
    fi
done