#!/bin/bash
#SBATCH --job-name=filter_10kb_array
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=1-16
#SBATCH --output=/scratch/au08019/reviopeanut/logs/filter_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/filter_%A_%a.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

# Runs each genotype as a separate array job (1–16)
module load seqtk/1.4-GCC-13.3.0
module load pigz/2.8-GCCcore-13.3.0

RAW_DIR=/scratch/au08019/reviopeanut
OUT_DIR=${RAW_DIR}/filteredfastq_min10kb
MINLEN=10000
mkdir -p "$OUT_DIR" /scratch/au08019/reviopeanut/logs

# --- Identify file for this array task ---
files=(${RAW_DIR}/*.fq.gz)
fq=${files[$SLURM_ARRAY_TASK_ID-1]}
base=$(basename "$fq")
out=${OUT_DIR}/${base%.fq.gz}.min${MINLEN}.fq.gz

echo "Job ID: $SLURM_ARRAY_JOB_ID  |  Task ID: $SLURM_ARRAY_TASK_ID"
echo "Filtering file: $fq"
echo "Output file:    $out"
echo "Keeping reads ≥ ${MINLEN} bp"

# --- Run filtering ---
seqtk seq -L $MINLEN "$fq" | pigz -p 4 > "$out"

echo "Done: $base"
echo "Completed at: $(date)"