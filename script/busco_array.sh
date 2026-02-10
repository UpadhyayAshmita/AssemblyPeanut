#!/bin/bash
#SBATCH --job-name=busco_array
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=36:00:00
#SBATCH --array=1-16
#SBATCH --output=/scratch/au08019/reviopeanut/logs/busco_%A_%a.out
#SBATCH --error=/scratch/au08019/reviopeanut/logs/busco_%A_%a.err

set -euo pipefail

# ---- project paths ----
WD="/scratch/au08019/reviopeanut"
GENOMEDIR="${WD}/scaffolded_pteran"

LIST="${GENOMEDIR}/genomes.list"
OUT="${WD}/busco_out"

mkdir -p "${WD}/logs" "${OUT}"

# ---- select genome for this array task ----
GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")
GENOME_PATH="${GENOMEDIR}/${GENOME}"

# sample name (remove .fa.gz)
SAMPLE=$(basename "${GENOME}" .fa.gz)

# ---- load BUSCO module ----
module purge
module load BUSCO/5.8.3-foss-2023a

# ---- BUSCO settings ----
LINEAGE="embryophyta_odb10"

# ---- temporary work directory ----
TMP="/lscratch/$USER/busco_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMP}"

FASTA="${TMP}/${SAMPLE}.fa"
gzip -dc "${GENOME_PATH}" > "${FASTA}"

# ---- run BUSCO ----
busco \
  -i "${FASTA}" \
  -o "${SAMPLE}" \
  -m genome \
  -l "${LINEAGE}" \
  -c "${SLURM_CPUS_PER_TASK}" \
  --out_path "${OUT}"

# ---- cleanup ----
rm -rf "${TMP}"
