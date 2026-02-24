#!/bin/bash
#SBATCH --job-name=map100kb
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.err

set -euo pipefail

module load minimap2/2.29-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

BASE=/scratch/au08019/reviopeanut/tetrasomy
REF=$BASE/ref/diploids_combined.fa

SAMPLE=$1
READS=$2

BAM=$BASE/map/${SAMPLE}.toDiploids.bam
WIN=$BASE/windows/diploids_combined.100kb.bed
OUT=$BASE/cov/${SAMPLE}.diploids_combined.100kb.meanDepth.tsv

mkdir -p $BASE/map $BASE/windows $BASE/cov $BASE/logs

echo "Mapping $SAMPLE ..."

# ---- MAP ----
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi "$REF" "$READS" \
  | samtools view -b -q 20 -F 0x100 -F 0x800 - \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "$BAM" -

samtools index "$BAM"

echo "Mapping done."

# ---- CREATE 100kb WINDOWS (only once) ----
if [ ! -f "$WIN" ]; then
    echo "Creating 100 kb windows..."
    samtools faidx "$REF"
    cut -f1,2 ${REF}.fai > $BASE/windows/diploids_combined.genome
    bedtools makewindows \
        -g $BASE/windows/diploids_combined.genome \
        -w 100000 \
        > "$WIN"
fi

# ---- DEPTH PER WINDOW ----
echo "Computing mean depth..."

bedtools coverage \
    -a "$WIN" \
    -b "$BAM" \
    -mean \
    > "$OUT"

echo "Done."
wc -l "$OUT"
head "$OUT"