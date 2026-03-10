#!/bin/bash

echo "=============================================="
echo " Phase 1: Compiling Individual PSD Databases  "
echo "=============================================="

# Recursively find all .sdf.gz files in all subdirectories (AA, AB, etc.)
find . -type f -name "*.sdf.gz" | while IFS= read -r input_file; do
    
    # Define the output filename by swapping the extension
    output_file="${input_file%.sdf.gz}.psd"
    
    # Check if the compiled .psd file already exists to avoid duplicate work
    if [[ -f "$output_file" ]]; then
        echo "[SKIPPED] Already compiled: $output_file"
    else
        echo "[BUILDING] Compiling $input_file ..."
        # Run the CDPKit creation tool
        psdcreate -i "$input_file" -o "$output_file"
    fi

done

echo ""
echo "=============================================="
echo " Phase 2: Merging into Master Database        "
echo "=============================================="

# Create a text file listing the exact paths of all the new .psd files
echo "Gathering list of compiled databases..."
find . -type f -name "*.psd" > list_of_databases.txt

# Merge them into one final, highly compressed binary index
echo "Fusing databases... This will utilize heavy RAM."
psdmerge -i list_of_databases.txt -o ZINC_InStock_Master.psd

echo ""
echo "=============================================="
echo " BUILD COMPLETE! Master File: ZINC_InStock_Master.psd "
echo "=============================================="