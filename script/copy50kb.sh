#!/bin/bash
#SBATCH --job-name=copy_zip_50k
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/logs/copy_zip_50k_%A.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/copy_zip_50k_%A.err


BASE_DIR=/scratch/au08019/reviopeanut/assembly
OUT_DIR=${BASE_DIR}/50k_files

mkdir -p "$OUT_DIR"

module load pigz/2.8-GCCcore-13.3.0   # parallel gzip 

echo "Starting copy and compression process..."

for dir in ${BASE_DIR}/*/; do
    cd "$dir" || continue
    GENO=$(basename "$dir")
    FILES=(*50k*.fa)

    # Check if any matching file exists
    if ls *50k*.fa >/dev/null 2>&1; then
        for f in *50k*.fa; do
            echo "Copying and compressing $GENO → $f"
            cp "$f" "${OUT_DIR}/${GENO}_$(basename "$f")"
            pigz -p 2 "${OUT_DIR}/${GENO}_$(basename "$f")"
        done
    else
        echo "No 50k .fa file found in $GENO"
    fi
done

echo "All matching files copied and compressed in: $OUT_DIR"
