#!/bin/bash
input_file="$1"
output_file="${input_file%.sdf.gz}.psd"

# Only compile if the output file doesn't already exist
if [[ ! -f "$output_file" ]]; then
    # Injecting the native thread limit directly into the software call
    psdcreate -i "$input_file" -o "$output_file" --num-threads 1 > /dev/null 2>&1
    echo "[DONE] $output_file"
fi