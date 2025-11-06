#!/bin/bash
#SBATCH --job-name=genomescope_peanut
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256gb
#SBATCH --time=36:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/output/genomescope/genomescope.%j.log
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL

module purge
module load Jellyfish/2.3.1-GCC-13.3.0
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

    #Count k-mers
    jellyfish count -C -m $k -s 56G -t 32 -o $sample_out/${sample}_${k}mer.jf <(zcat $fq)

    # Make histogram
    jellyfish histo -t 32 $sample_out/${sample}_${k}mer.jf > $sample_out/${sample}_${k}mer.histo

    #Run GenomeScope2 (the module provides the R package + CLI)
    genomescope2 \
        -i $sample_out/${sample}_${k}mer.histo \
        -k $k \
        -o $sample_out \
        -g $GENOME_SIZE
done