#!/bin/bash

# Define your input query and your final output file
QUERY="/mnt/e/TFEB_hopping/pharmacophore_west_face1.pml"
MASTER_HITS="All_Fragment_Hits.sdf"

# 1. Clean up temporary files from any previous runs
rm -f temp_hits.sdf "$MASTER_HITS"

echo "Initiating Fragmented Screening Architecture..."

# 2. Find all individual tranche databases, explicitly ignoring the massive master file
find . -type f -name "*.psd" ! -name "ZINC_InStock_Master.psd" | while read -r db_file; do
    
    # 3. Screen the tiny file (Loads instantly, blasts 46 cores, flushes RAM)
    psdscreen -Q pml -q "$QUERY" -i "$db_file" -o temp_hits.sdf -x 1 > /dev/null 2>&1
    
    # 4. If hits were found (the temp file is not empty), safely append them
    if [[ -s temp_hits.sdf ]]; then
        cat temp_hits.sdf >> "$MASTER_HITS"
        echo "[HITS FOUND] Appended matches from $db_file"
    fi
    
done

# 5. Clean up the temp file
rm -f temp_hits.sdf

echo "=========================================================="
echo "Screening Complete! Master hit list saved to: $MASTER_HITS"
echo "=========================================================="