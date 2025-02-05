#!/bin/bash
#SBATCH --partition=bahl_p
#SBATCH --job-name=BEAST
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-user=ll22780@uga.edu
#SBATCH --mail-type=END,FAIL

# Load Beast module
module load Beast/1.10.4-GCC-11.3.0

# Set the base directory
base_dir="./output"

# Create a directory named 'combined' under $base_dir
mkdir -p "${base_dir}/combined"

# Get the list of run directories
run_dirs=($(ls -d ${base_dir}/run*/))

# Get the list of log files in the first run directory
log_files=($(find "${run_dirs[0]}" -name "*.log"))

# Iterate over each log file in the first run directory
for log_file in "${log_files[@]}"; do
  file_name=$(basename "$log_file")
  found=true
  paths=()
  command="logcombiner -burnin 100000000 -resample 1600000"

  # Check if the log file exists in the other run directories
  for run_dir in "${run_dirs[@]}"; do
    if [ -f "${run_dir}${file_name}" ]; then
      paths+=("${run_dir}${file_name}")
    else
      found=false
      break
    fi
  done

  # If the log file was found in all run directories, print its paths in one line
  if [ "$found" = true ]; then
    combined_file="${base_dir}/combined/combined_$file_name"
    command+=" ${paths[@]} $combined_file"
    eval "$command"
  fi

done

# Get the list of trees files in the first run directory
trees_files=($(find "${run_dirs[0]}" -name "*.trees"))

# Iterate over each log file in the first run directory
for t_file in "${trees_files[@]}"; do
  file_name=$(basename "$t_file")
  found=true
  paths=()
  command="logcombiner -trees -burnin 100000000 -resample 1600000"

  # Check if the log file exists in the other run directories
  for run_dir in "${run_dirs[@]}"; do
    if [ -f "${run_dir}${file_name}" ]; then
      paths+=("${run_dir}${file_name}")
    else
      found=false
      break
    fi
  done

  # If the log file was found in all run directories, print its paths in one line
  if [ "$found" = true ]; then
    combined_file="${base_dir}/combined/combined_$file_name"
    command+=" ${paths[@]} $combined_file"
    eval "$command"
  fi

done
