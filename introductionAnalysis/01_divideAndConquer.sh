#!/bin/bash

# Fixed file paths
input_file="data/combined.trees"
output_dir="tree"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Initialize variables
header=""
footer=""
in_tree_block=0
tree_count=0

# Use awk to process the file in a single pass
awk -v header_ref="$header" -v footer_ref="$footer" -v output_dir="$output_dir" '
BEGIN {
    header = header_ref;
    footer = footer_ref;
    in_tree_block = 0;
    tree_count = 0;
}
{
    if ($0 ~ /^tree STATE/) {
        if (in_tree_block == 1) {
            # Write the previous tree content
            tree_count++;
            output_file = output_dir "/tree" tree_count ".nexus";
            print header tree_content footer "End;" > output_file;
            tree_content = "";
        }
        in_tree_block = 1;
    }
    if (in_tree_block == 0) {
        header = header $0 "\n";
    } else {
        tree_content = tree_content $0 "\n";
    }
}
END {
    if (in_tree_block == 1) {
        # Write the last tree content
        tree_count++;
        output_file = output_dir "/tree" tree_count ".nexus";
        print header tree_content footer > output_file;
    }
    print "Total trees extracted: " tree_count;
}' "$input_file"


# Loop through each generated Nexus file to create and submit a job script
for file in "$output_dir"/*.nexus; do
    tree="$(basename "$file" .nexus)"
    job_script="$output_dir"/"$tree"ToTB.sh

    # Create the job script
    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --partition=bahl_p
#SBATCH --job-name=RtransferTreeToTable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=48:00:00
#SBATCH --mail-user=ll22780@uga.edu
#SBATCH --mail-type=FAIL
module load R/4.3.1-foss-2022a

Rscript ../scripts/importTreeSummary.R $(basename "$file") Texas 2021.81369863014
EOT
    cd "$output_dir"
    sbatch $(basename "$job_script")
    cd ..
done

# Final message
echo "All jobs submitted"
