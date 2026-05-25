#!/bin/bash

src_dirs=(
    "/lustre/home/harrell_lab/WGS/processed/SUS20240514123-WGS1/BQSR"
    "/lustre/home/harrell_lab/WGS/processed/SUS20241016061-WGS2/BQSR"
)
target_dir="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS"

mkdir -p "$target_dir"

for src in "${src_dirs[@]}"; do
    find -L "$src" -type f | while read -r file; do
        # Resolve the real path to get the physical location
        real_src=$(realpath "$file")
        
        # Compute relative path from the original base dir
        rel_path="${file#$src/}"

        # Final symlink destination
        dest="$target_dir/$rel_path"

        # Ensure the target subdirectory exists
        mkdir -p "$(dirname "$dest")"

        # Create the symlink
        ln -sf "$real_src" "$dest"
    done
done
