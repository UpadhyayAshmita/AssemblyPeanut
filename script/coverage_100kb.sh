#!/bin/bash
#SBATCH --job-name=cov100kb_TifTB
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.err

set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0

cd /scratch/au08019/reviopeanut/tetrasomy

# inputs
BED=windows/diploids_combined.100kb.bed
BAM=map/TifTB.toDiploids.bam

# output
OUT=cov/TifTB.diploids_combined.100kb.meanDepth.tsv

mkdir -p cov logs

echo "[$(date)] Starting bedtools coverage"
echo "BED: $BED"
echo "BAM: $BAM"
echo "OUT: $OUT"

# sanity checks
test -s "$BED"
test -s "$BAM"
test -s "${BAM}.bai" || samtools index "$BAM"

bedtools coverage \
  -a "$BED" \
  -b "$BAM" \
  -mean \
  > "$OUT"

echo "[$(date)] Done"
wc -l "$OUT"
head "$OUT"