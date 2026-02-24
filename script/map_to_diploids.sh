#!/bin/bash
#SBATCH --job-name=window_depths
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/tetrasomy/logs/%x.%j.err

set -euo pipefail

# ---- modules ----
module load minimap2/2.29-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0

# ---- inputs ----
BASE=/scratch/au08019/reviopeanut/tetrasomy
REF=$BASE/ref/diploids_chr10_combined.fa

SAMPLE=TifTB
READS=/scratch/au08019/reviopeanut/TifTB.fq.gz

# ---- outputs ----
mkdir -p $BASE/map $BASE/logs
BAM=$BASE/map/${SAMPLE}.toDiploids.bam

echo "[$(date)] Mapping $SAMPLE to diploid ref..."
echo "REF:   $REF"
echo "READS: $READS"
echo "OUT:   $BAM"

# # ---- map -> sort ----
# minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi "$REF" "$READS" \
#   | samtools view -b -q 20 -F 0x100 -F 0x800 - \
#   | samtools sort -@ 16 -o "$BAM" -

# # ---- index + stats ----
# samtools index "$BAM"
# samtools flagstat "$BAM" > $BASE/logs/${SAMPLE}.toDiploids.flagstat.txt

# echo "[$(date)] Done."
# tail -n 12 $BASE/logs/${SAMPLE}.toDiploids.flagstat.txt


module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0

BASE=/scratch/au08019/reviopeanut/tetrasomy

# echo "Starting coverage calculation at $(date)"
mkdir -p $BASE/windows

cut -f1,2 $REF.fai > $BASE/windows/diploids_allchr.genome

BASE=/scratch/au08019/reviopeanut/tetrasomy
REF=$BASE/ref/diploids_combined.fa

mkdir -p $BASE/windows

# genome file (contig  length)
cut -f1,2 ${REF}.fai > $BASE/windows/diploids_combined.genome

# 100 kb windows
bedtools makewindows \
  -g $BASE/windows/diploids_combined.genome \
  -w 100000 \
  > $BASE/windows/diploids_combined.100kb.bed

wc -l $BASE/windows/diploids_combined.100kb.bed
head $BASE/windows/diploids_combined.100kb.bed