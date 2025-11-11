#!/bin/bash
#SBATCH --job-name=quast_primary_ref
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=36:00:00
#SBATCH --array=1-16
#SBATCH --output=/scratch/au08019/reviopeanut/logs/quast_ref_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/quast_ref_%A_%a.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

module purge
module load QUAST/5.2.0-gfbf-2023b

ASM_DIR=/scratch/au08019/reviopeanut/assembly
OUT_DIR=/scratch/au08019/reviopeanut/quast_primary_ref
REF=/scratch/au08019/reviopeanut/reference/tifrunner_v2.fa

mkdir -p $OUT_DIR
folders=(${ASM_DIR}/*)
FOLDER=${folders[$SLURM_ARRAY_TASK_ID-1]}
BASENAME=$(basename $FOLDER)

ASM_FILE=${FOLDER}/${BASENAME}.asm.bp.p_ctg.fa
OUT_SUB=${OUT_DIR}/${BASENAME}
mkdir -p $OUT_SUB
echo "==============================================================="
echo "Running QUAST (reference-based) for $BASENAME"
echo "Input assembly: $ASM_FILE"
echo "Reference: $REF"
echo "Output folder: $OUT_SUB"
echo "==============================================================="

if [ -f "$ASM_FILE" ]; then
    quast.py "$ASM_FILE" \
        -o "$OUT_SUB" \
        -r "$REF" \
        -t 16 
else
    echo "Assembly not found: $ASM_FILE"
fi

echo "Done with $BASENAME"
