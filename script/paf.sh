#!/bin/bash
#SBATCH --job-name=makePAF
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=04:00:00
#SBATCH --array=1-4
#SBATCH --output=/scratch/au08019/reviopeanut/logs/makePAF_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/makePAF_%A_%a.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

module load minimap2/2.28-GCCcore-13.2.0

REF="/scratch/au08019/reviopeanut/reference/tifrunner_v2_filt.fa"
IN="/scratch/au08019/reviopeanut/assembly/broken_genome"
OUT="/scratch/au08019/reviopeanut/assembly/paf"
GENOFILE="/scratch/au08019/reviopeanut/assembly/broken_genome/genos.txt"

mkdir -p "$OUT"

# get Nth line
FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENOFILE")
BASE="${FA%.fa.gz}"; BASE="${BASE%.fa}"

echo "Running $FA"

minimap2 -x asm5 -t 32 "$REF" "$IN/$FA" > "$OUT/${BASE}_vs_Tifrunner.paf"

echo "Done: $OUT/${BASE}_vs_Tifrunner.paf"
