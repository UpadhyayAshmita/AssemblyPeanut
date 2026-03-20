#!/bin/bash
#SBATCH --job-name=mm2_broken3
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --output=/scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/08_logs/mm2_broken3_%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/08_logs/mm2_broken3_%j.err

set -euo pipefail

echo "=== Job started ==="
date
hostname

module purge
module load minimap2

echo "Minimap2 path:"
which minimap2

REF="/scratch/au08019/reviopeanut/reference/tifrunner_v2_filt.fa"
QUERY="/scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/06_fasta_fix/output/TifNV-HG.asm.bp.p_ctg.broken3.fa"
OUT="/scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/07_plots/TifNV-HG.asm.bp.p_ctg.broken3_vs_Tifrunner.paf"

echo "Reference: $REF"
echo "Query: $QUERY"
echo "Output: $OUT"

[ -f "$REF" ] || { echo "ERROR: Reference not found: $REF"; exit 1; }
[ -f "$QUERY" ] || { echo "ERROR: Query not found: $QUERY"; exit 1; }

mkdir -p /scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/07_plots
mkdir -p /scratch/au08019/reviopeanut/assembly/paf/TifNV_audit/08_logs

minimap2 -x asm5 --eqx -t "${SLURM_CPUS_PER_TASK}" \
  "$REF" \
  "$QUERY" \
  > "$OUT"

echo "=== Alignment done ==="
date
ls -lh "$OUT"