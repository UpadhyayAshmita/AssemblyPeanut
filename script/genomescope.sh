#!/bin/bash
#SBATCH --job-name=genomescope_peanut
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/output/genomescope/genomescope.%j.log
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

module load Jellyfish/2.3.0-GCC-12.3.0
module load R/4.3.2
module load GenomeScope/2.0.1-foss-2023a-R-4.3.2

READS=/scratch/au08019/reviopeanut
OUT=/scratch/au08019/reviopeanut/output/genomescope
mkdir -p $OUT

k=21
GENOME_SIZE=2.7e9

echo "Job started at: $(date)"
for fq in $READS/*.fq.gz; do
    sample=$(basename ${fq%.fq.gz})
    echo "[$(date)] Processing $sample ..."

    sample_out=$OUT/${sample}
    mkdir -p $sample_out

    # Step 1: Count k-mers
    jellyfish count -C -m $k -s 100G -t 32 -o $sample_out/${sample}_${k}mer.jf <(zcat $fq)

    # Step 2: Create histogram
    jellyfish histo -t 32 $sample_out/${sample}_${k}mer.jf > $sample_out/${sample}_${k}mer.histo

    # Step 3: Run GenomeScope2
    Rscript -e "library(genomescope2.0); genomescope2.0('${sample_out}/${sample}_${k}mer.histo', $k, '${sample_out}', $GENOME_SIZE)"
done

echo "[$(date)] GenomeScope2 analysis complete. Results in: $OUT"
