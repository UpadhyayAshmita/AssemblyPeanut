#!/bin/bash
#SBATCH --job-name=hifiasm_array
#SBATCH --partition=batch
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --mem=512G
#SBATCH --time=72:00:00
#SBATCH --array=1-16
#SBATCH --output=/scratch/au08019/reviopeanut/logs/hifiasm_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/hifiasm_%A_%a.err

module load hifiasm/0.25.0

WORKDIR=/scratch/au08019/reviopeanut
READLIST=${WORKDIR}/sample_fastqs.txt

# Get the fastq for this array index
FQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$READLIST")
BASENAME=$(basename "$FQ" .min10000.fq.gz)

OUTDIR=${WORKDIR}/assembly/${BASENAME}
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Assembling sample: $BASENAME"
echo "Reads: $FQ"

hifiasm -t ${SLURM_CPUS_PER_TASK} \
  -o ${BASENAME}.asm \
  "$FQ" 2> ${BASENAME}.asm.log

# Convert primary contigs to FASTA
awk '/^S/{print ">"$2"\n"$3}' ${BASENAME}.asm.p_ctg.gfa > ${BASENAME}.asm.p_ctg.fa
