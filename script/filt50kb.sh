#!/bin/bash
#SBATCH --job-name=filter50k
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-16
#SBATCH --output=/scratch/au08019/reviopeanut/logs/filter50k_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/filter50k_%A_%a.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

# Filter contigs <50 kb for all assemblies (HiFiasm primary)

module purge
module load seqtk/1.4-GCC-13.3.0

ASM_DIR=/scratch/au08019/reviopeanut/assembly

# Collect list of assembly folders
folders=(${ASM_DIR}/*)
FOLDER=${folders[$SLURM_ARRAY_TASK_ID-1]}
BASENAME=$(basename $FOLDER)

INPUT=${FOLDER}/${BASENAME}.asm.bp.p_ctg.fa
OUTPUT=${FOLDER}/${BASENAME}.asm.bp.p_ctg.50k.fa

echo "Filtering contigs >=50kb for: $BASENAME"
echo "Input:  $INPUT"
echo "Output: $OUTPUT"

if [ -f "$INPUT" ]; then
    seqtk seq -L 50000 "$INPUT" > "$OUTPUT"

    echo "Counting sequences..."
    echo "Before filtering: $(grep -c '^>' $INPUT)"
    echo "After filtering:  $(grep -c '^>' $OUTPUT)"
else
    echo " Missing input FASTA: $INPUT"
fi

echo "Done with $BASENAME"
