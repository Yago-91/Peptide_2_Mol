#!/bin/bash

# 1. Setup global variables
export QUERY="YOUR_MODEL_NAME.pml"
export FINAL_HITS="All_Fragment_Hits.sdf"

# 2. The Master Lock: Force each psdscreen instance to use exactly 1 thread
export OMP_NUM_THREADS=1

# Clean up any previous runs to start fresh
rm -f "$FINAL_HITS"
find . -type f -name "*_hits.sdf" -delete

echo "Building the worker script..."

# 3. Create the "Worker" script on the fly
cat << 'EOF' > worker_screen.sh
#!/bin/bash
db_file="$1"
# Create a unique output file name for EVERY database (e.g., AA_AARN_hits.sdf)
output_file="${db_file%.psd}_hits.sdf"

# Run the screen silently
psdscreen -Q pml -q "$QUERY" -d "$db_file" -o "$output_file" -x 1 -t 1 > /dev/null 2>&1

# If the output file is totally empty (no hits found), delete it to keep things clean
if [[ ! -s "$output_file" ]]; then
    rm -f "$output_file"
fi
EOF

chmod +x worker_screen.sh

echo "Deploying 46 independent screening jobs simultaneously..."

# 4. The Scatter: Feed the small databases to xargs, maintaining 46 active jobs
find . -type f -name "*.psd" ! -name "ZINC_InStock_Master.psd" | xargs -n 1 -P 46 ./worker_screen.sh

echo "Screening complete! Gathering and merging hit files..."

# 5. The Gather: Safely concatenate all the unique hit files into one master file
# Using find and exec bypasses the Linux argument limit safely
find . -type f -name "*_hits.sdf" -exec cat {} + > "$FINAL_HITS"

echo "=========================================================="
echo "Master hit list successfully generated: $FINAL_HITS"
echo "=========================================================="