#!/bin/bash
#SBATCH --job-name=gfa2fasta
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/logs/gfa2fasta_%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/gfa2fasta_%j.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL


WORKDIR=/scratch/au08019/reviopeanut/assembly

echo "Starting GFA → FASTA conversion..."
echo "Working directory: $WORKDIR"
echo "---------------------------------------------"

# Loop through all assembly folders
for d in ${WORKDIR}/*/; do
    cd "$d" || continue
    BASENAME=$(basename "$d")

    echo "Processing $BASENAME ..."

    # Find all .p_ctg.gfa, hap1, hap2, etc.
    for f in ${BASENAME}.asm*.p_ctg.gfa; do
        if [ -f "$f" ]; then
            OUT="${f%.gfa}.fa"
            echo "   Converting $f → $(basename $OUT)"
            awk '/^S/{print ">"$2"\n"$3}' "$f" > "$OUT"
        else
            echo "No file found for $f"
        fi
    done

    echo " Done: $BASENAME"
    echo "---------------------------------------------"
done

echo "All conversions completed!"
