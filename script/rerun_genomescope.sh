#!/bin/bash
#SBATCH --job-name=genomescope_rerun
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/output/genomescope/rerun_genomescope.%j.log
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL


module purge
module load R/4.3.2

export R_LIBS_USER="/home/au08019/R/x86_64-pc-linux-gnu-library/4.3"

# Force headless plotting
export R_DEFAULT_DEVICE="png"

HISTO_DIR="/scratch/au08019/reviopeanut/output/genomescope"
GENOMESCOPE_SCRIPT="/home/au08019/genomescope2.0/genomescope.R"

echo "Job started: $(date)"
cd "$HISTO_DIR"

# Loop through existing histograms (*.histo)

for histo in */*_21mer.histo; do
    sample=$(basename "$(dirname "$histo")")

    echo "[$(date)] Running GenomeScope for $sample"

    Rscript --vanilla "$GENOMESCOPE_SCRIPT" \
        -i "$HISTO_DIR/$sample/${sample}_21mer.histo" \
        -o "$HISTO_DIR/$sample" \
        -k 21

    echo "[$(date)] Completed GenomeScope for $sample"
done

echo "Job finished: $(date)"
