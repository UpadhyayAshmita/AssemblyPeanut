#!/bin/bash
#SBATCH --job-name=copy_zip_50k_missing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/logs/copy_zip_missing_%A.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/copy_zip_missing_%A.err

BASE_DIR=/scratch/au08019/reviopeanut/assembly
OUT_DIR=${BASE_DIR}/50k_files

mkdir -p "$OUT_DIR"

module load pigz/2.8-GCCcore-13.3.0

echo "Checking for missing 50k files..."

for dir in ${BASE_DIR}/*/; do
    GENO=$(basename "$dir")
    
    # Skip output folder itself
    if [[ "$GENO" == "50k_files" ]]; then
        continue
    fi

    cd "$dir" || continue

    # Check if the folder has any 50k FASTA
    FILES=( *50k*.fa )
    if [[ ! -e "${FILES[0]}" ]]; then
        echo "No 50k .fa file found in $GENO"
        continue
    fi

    for f in *50k*.fa; do
        OUT_FILE="${OUT_DIR}/${GENO}_$(basename "$f").gz"

        # Skip if gzipped output already exists
        if [[ -f "$OUT_FILE" ]]; then
            echo "Already done: $OUT_FILE (skipping)"
            continue
        fi

        echo "Copying and compressing $GENO → $f"
        cp "$f" "${OUT_DIR}/${GENO}_$(basename "$f")"
        pigz -p 2 "${OUT_DIR}/${GENO}_$(basename "$f")"
    done
done

echo "Done. Only missing files processed."
