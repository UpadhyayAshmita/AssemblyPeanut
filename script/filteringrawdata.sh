#!/bin/bash
#SBATCH --job-name=filter_hifi_10kb
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/au08019/reviopeanut/logs/filter_hifi_10kb.%j.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/filter_hifi_10kb.%j.err
#SBATCH --mail-user=au08019@uga.edu
#SBATCH --mail-type=END,FAIL
# Script: Filter PacBio HiFi reads ≥10 kb using seqtk


module load seqtk/1.4-GCC-13.3.0

RAW_DIR=/scratch/au08019/reviopeanut
OUT_DIR=${RAW_DIR}/filteredfastq_min10kb
MINLEN=10000   # 10 kb = 10,000 bp

mkdir -p "${OUT_DIR}"

echo "Filtering all .fq.gz files in ${RAW_DIR} to keep reads ≥${MINLEN} bp..."
echo "Output will be saved in: ${OUT_DIR}"
for fq in ${RAW_DIR}/*.fq.gz; do
    base=$(basename "${fq}")
    out="${OUT_DIR}/${base%.fq.gz}.min${MINLEN}.fq.gz"
    echo "Processing: ${base}"
    seqtk seq -L ${MINLEN} "${fq}" | gzip > "${out}"
    echo "Done -> ${out}"
done
echo "All files filtered successfully."
