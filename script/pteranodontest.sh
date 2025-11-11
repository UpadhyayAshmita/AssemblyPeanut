#!/bin/bash
#SBATCH --job-name=ptera_15X034_50k
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/logs/ptera_15X034_50k_%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/ptera_15X034_50k_%j.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

#  Pteranodon scaffolding for one assembly (50 kb filtered)

module purge
module load GCC/13.2.0
module load Miniconda3/24.1.2-0
source activate pteranodon_env

# --- Input / Output paths ---
REF=/scratch/au08019/reviopeanut/reference/tifrunner_v2.fa
ASM=/scratch/au08019/reviopeanut/assembly/15X034-1-1-SSD-15/15X034-1-1-SSD-15.asm.bp.p_ctg.50k.fa
OUTDIR=/scratch/au08019/reviopeanut/scaffolded/15X034-1-1-SSD-15

mkdir -p $OUTDIR

echo "Running Pteranodon scaffolding"
echo "Reference : $REF"
echo "Assembly  : $ASM"
echo "OutputDir : $OUTDIR"

bash ~/apps/Pteranodon/PteranodonBase.sh \
    -ref $REF \
    -query $ASM \
    -o $OUTDIR/15X034-1-1-SSD-15 \
    -t 32 \
    -SegLen 20000 \
    -MinQueryLen 0.05 \
    --auto

echo "Scaffolding completed for 15X034-1-1-SSD-15 (50 kb filtered)"
