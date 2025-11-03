#!/bin/bash
#SBATCH --job-name=hifi_QC
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/output/qc_reports/hifi_QC.%j.log
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL


module load FastQC/0.11.9-Java-11
module load MultiQC/1.28-foss-2024a
module load SeqKit/2.9.0
module load NanoPlot/1.43.0-foss-2023a

Reads=/scratch/au08019/reviopeanut
OUT=/scratch/au08019/reviopeanut/output/qc_reports

mkdir -p $OUT/fastqc $OUT/nanoplot $OUT/seqkit

# Run FastQC on all HiFi reads
echo "Running FastQC..."
fastqc -t 16 -o $OUT/fastqc $Reads/*.fq*


# Generate NanoPlot QC for read length & quality

echo "Running NanoPlot..."
for f in $Reads/*.fq*; do
    base=$(basename ${f%.*})
    NanoPlot --fastq $f \
             -t 8 \
             --maxlength 60000 \
             --plots hex dot \
             --outdir $OUT/nanoplot/${base} \
             --title "${base}_NanoPlot"
done


# Run SeqKit stats for each file
echo "Generating SeqKit statistics..."
for f in $Reads/*.fq*; do
    base=$(basename ${f%.*})
    seqkit stats $f > $OUT/seqkit/${base}_stats.txt
done
# Summarize all results using MultiQC
echo "Aggregating QC reports with MultiQC..."
multiqc $OUT -o $OUT
echo "✅ QC complete! All reports are in: $OUT"
