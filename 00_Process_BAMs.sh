# Load settings
source ./00_Settings.sh

# Create necessary directories
mkdir -p ./Processed_BAMs
mkdir -p ./sorted

# Loop through BAM files, following symbolic links
find -L "$BAM_DIR" -type f -name "$BAM_REGEX" | while read -r BAM; do
    BAM_base=$(basename "$BAM")
    BAM_sorted="./sorted/$(basename "$BAM" .bam).sorted.bam"

    # Sort and index the BAM file
    samtools sort -@ "$MAX_CPUS" "$BAM" -o "$BAM_sorted"
    samtools index -@ "$MAX_CPUS" "$BAM_sorted"

    # Create a list of valid chromosomes (those that appear in human_rename.txt)
    valid_chrs=()
    while read -r line; do
        valid_chrs+=("$(echo "$line" | awk '{print $1}')")
    done < human_rename.txt

    # Filter the BAM file for valid chromosomes
    samtools view -@ "$MAX_CPUS" -b "$BAM_sorted" "${valid_chrs[@]}" > temp_filtered.bam

    # Extract and modify the header (rename chromosomes in the header)
    samtools view -H temp_filtered.bam | \
    awk -v OFS="\t" '
        BEGIN { 
            # Build the mapping from old name to new name from "human_rename.txt"
            while ((getline < "human_rename.txt") > 0) 
                a["SN:" $1] = "SN:" $2 
        }
        /^@SQ/ { 
            # Loop through each field to find the "SN:" field
            for (i=1; i<=NF; i++) 
                if ($i ~ /^SN:/) {
                    # If the sequence name matches the mapping, replace it
                    if ($i in a) {
                        sub($i, a[$i]); 
                        print $0; 
                        next;
                    }
                }
            # Skip printing @SQ lines that do not match any in the mapping
            next;
        }
        # Print all non-@SQ lines as is
        { print }
    ' > new_header.sam

    # Replace the header in the filtered BAM file
    samtools reheader new_header.sam temp_filtered.bam > temp_reheadered.bam
    samtools sort -@ "$MAX_CPUS" temp_reheadered.bam -o "./Processed_BAMs/${BAM_base}"

    samtools index -@ "$MAX_CPUS" "./Processed_BAMs/${BAM_base}"

    # Clean up temporary files
    rm temp_filtered.bam temp_reheadered.bam new_header.sam
done
