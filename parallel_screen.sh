#!/bin/bash

# Default to 40 threads if the user doesn't specify one
THREADS=40

# Help menu function
usage() {
    echo "Usage: $0 -q <query.pml> -o <output.sdf> [-t <threads>]"
    echo "  -q : Path to your LigandScout .pml file (Required)"
    echo "  -o : Name of the final output .sdf file (Required)"
    echo "  -t : Number of threads/concurrent jobs to run (Optional, default 40)"
    exit 1
}

# Parse command-line arguments
while getopts "q:o:t:h" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        o) FINAL_HITS="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Ensure the required arguments were provided
if [ -z "$QUERY" ] || [ -z "$FINAL_HITS" ]; then
    echo "Error: Missing required arguments."
    usage
fi

if [ ! -f "$QUERY" ]; then
    echo "Error: Cannot find pharmacophore file '$QUERY'."
    exit 1
fi

# Create a unique Run ID using the script's Process ID ($$) to prevent collisions
# if you decide to run multiple different screens at the same time.
RUN_ID="run_$$"
LOG_FILE="${RUN_ID}_completed.log"
WORKER_SCRIPT="${RUN_ID}_worker.sh"

rm -f "$FINAL_HITS" "$LOG_FILE"
touch "$LOG_FILE"

echo "Scanning directory to count databases..."
TOTAL_FILES=$(find . -type f -name "*.psd" ! -name "ZINC_InStock_Master.psd" | wc -l)
echo "Found $TOTAL_FILES databases. Initializing workload..."

# Build the worker script (using EOF without quotes so we can pass the RUN_ID and QUERY)
cat << EOF > "$WORKER_SCRIPT"
#!/bin/bash
db_file="\$1"
# Create a completely unique temp file name for this specific run
output_file="./${RUN_ID}_temp_\$(basename "\$db_file" .psd).sdf"

psdscreen -Q pml -q "$QUERY" -d "\$db_file" -o "\$output_file" -x 1 -t 1 > /dev/null 2>&1

if [[ ! -s "\$output_file" ]]; then
    rm -f "\$output_file"
fi

echo "." >> "$LOG_FILE"
EOF

chmod +x "$WORKER_SCRIPT"

# The Background Monitor Function
monitor_progress() {
    while true; do
        DONE=\$(wc -l < "$LOG_FILE")
        if [ "$TOTAL_FILES" -gt 0 ]; then
            PERCENT=\$(( DONE * 100 / TOTAL_FILES ))
        else
            PERCENT=100
        fi
        
        printf "\rProgress: [%-50s] %d%% (%d/%d)" \$(printf "#%.0s" \$(seq 1 \$((PERCENT / 2)))) "\$PERCENT" "\$DONE" "\$TOTAL_FILES"
        
        if [ "\$DONE" -ge "$TOTAL_FILES" ]; then break; fi
        sleep 2
    done
}

echo "Deploying $THREADS independent screening threads for $QUERY..."

monitor_progress &
MONITOR_PID=$!

# The Scatter
find . -type f -name "*.psd" ! -name "ZINC_InStock_Master.psd" | xargs -n 1 -P "$THREADS" ./"$WORKER_SCRIPT"

wait $MONITOR_PID
echo "" 

echo "Gathering and merging hit files..."

# The Gather: Only collect temp files belonging to THIS specific run
find . -maxdepth 1 -type f -name "${RUN_ID}_temp_*.sdf" -exec cat {} + > "$FINAL_HITS"

# Clean up ALL temp files and trackers associated with this run
rm -f "$LOG_FILE" "$WORKER_SCRIPT"
find . -maxdepth 1 -type f -name "${RUN_ID}_temp_*.sdf" -delete

echo "=========================================================="
echo "Master hit list successfully generated: $FINAL_HITS"
echo "=========================================================="