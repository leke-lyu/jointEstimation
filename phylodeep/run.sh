#!/bin/bash

# Directory containing the .tre files
TREE_DIR="./tre"

# Activate conda environment
source ~/.bashrc
conda activate phyloenv

# Loop through all .tre files
for filepath in "$TREE_DIR"/*.tre; do
    echo "Processing $filepath"

    # Step 1: Extract the Newick string from the tree line
    newick=$(grep -E '^\s*tree' "$filepath" | sed -E 's/.*=\s*\[&R\]\s*//')

    # Check that we got something before overwriting
    if [ -n "$newick" ]; then
        echo "$newick" > "$filepath"
    else
        echo "Warning: No tree line found in $filepath"
    fi

    # Step 2: Run paramdeep
    output="${filepath%.tre}_BD_FFNN.csv"
    paramdeep -t "$filepath" -p 0.1 -m BD -v FFNN_SUMSTATS -o "$output"
done

echo "Done processing all .tre files."
