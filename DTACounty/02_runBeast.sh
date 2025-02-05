files=$(find . -maxdepth 1 -type f -name "*.xml")
  
for file in $files
do
  base=$(basename "$file" .xml)
  if [ ! -d "$base" ]; then # Create a directory based on the base name if it doesn't already exist
    mkdir "$base"
    echo "Directory created: $base"
    cd "$base"
        for i in $(seq 1 5); do
                mkdir "run$i"
                cd run$i
                cat << EOF > "run.sh"
#!/bin/bash
#SBATCH --partition=bahl_p
#SBATCH --nodelist=ra1-18
#SBATCH --job-name=BEAST
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2gb
#SBATCH --time=720:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-user=ll22780@uga.edu
#SBATCH --mail-type=END,FAIL
module load Beast/1.10.4-GCC-11.3.0
beast ../../$file
EOF
                sbatch run.sh
                cd ..
        done
    cd ..
  else
    echo "Directory already exists: $base"
  fi
done
