#!/bin/bash

# Define the directory containing the .trees files
input_dir="output/combined"

# Iterate over all .trees files in the specified directory
for input_file in "$input_dir"/*.trees; do
  # Get the base name of the file (without the .trees extension)
  base_name=$(basename "$input_file" .trees)
  
  # Initialize variables
  header=""
  footer=""
  in_tree_block=0
  tree_count=0

  echo "Processing file: $input_file with base name: $base_name"

  # Use awk to process the file in a single pass
  awk -v header_ref="$header" -v footer_ref="$footer" -v base_name="$base_name" -v output_dir="$input_dir" '
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
              output_file = output_dir "/" tree_count "_" base_name ".tree";
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
          output_file = output_dir "/" tree_count "_" base_name ".tree";
          print header tree_content footer > output_file;
      }
      print "Total trees extracted from " base_name ": " tree_count > "/dev/stderr";  # Debugging statement
  }' "$input_file"
done

cd "$input_dir"
# Iterate over sequence numbers from 1 to 500
for i in $(seq 1 1000); do
  # Find all files that start with i_
  files=$(ls ${i}_*.tree 2>/dev/null)
  
  # Check if any files were found
  if [ -n "$files" ]; then
    # Combine the files into a single input for the R script
    combined_files=""
    for file in $files; do
      combined_files="$combined_files $file"
    done
    
    job_script="$i".sh
    # Create the job script
    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --partition=bahl_p
#SBATCH --job-name=sss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=48:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-user=ll22780@uga.edu
#SBATCH --mail-type=FAIL
module load R/4.3.1-foss-2022a

Rscript ../../scripts/processTrees.R ${combined_files} ${i}.RData
EOT

   sbatch "$job_script"

  fi
done
