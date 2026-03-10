#!/bin/bash
input_file="$1"
output_file="${input_file%.sdf.gz}.psd"

# Only compile if the output file doesn't already exist
if [[ ! -f "$output_file" ]]; then
    # Silently run the compilation to keep the terminal clean
    psdcreate -i "$input_file" -o "$output_file" > /dev/null 2>&1
    echo "[DONE] $output_file"
fi