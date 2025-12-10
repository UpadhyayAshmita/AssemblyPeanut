#!/bin/bash
#SBATCH --job-name=hifi_QC
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/output/qc_reports/hifi_QC.%j.log
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL


module load FastQC/0.11.9-Java-11
module load MultiQC/1.28-foss-2024a
module load NanoPlot/1.43.0-foss-2023a

Reads=/scratch/au08019/reviopeanut
OUT=/scratch/au08019/reviopeanut/output/qc_reports

mkdir -p $OUT/nanoplot 

# Generate NanoPlot QC for read length & quality
f="$Reads/C2563-7.fq.gz"  
base=$(basename ${f%.*})

echo "Running NanoPlot for: $base"

NanoPlot --fastq $f \
         -t 32 \
         --maxlength 60000 \
         --plots hex dot \
         --outdir $OUT/nanoplot/${base} \
         --title "${base}_NanoPlot"
