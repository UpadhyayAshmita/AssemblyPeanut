#!/bin/bash
#SBATCH -J ptero_single
#SBATCH --time=48:00:00
#SBATCH -c 24
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=180G
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%j.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%j.err"

# Load Khufu + Pteranodon env
source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnv.sh
source /cluster/projects/khufu/korani_projects/load_modules.sh

BASE_DIR="/cluster/lab/clevenger/Ashmita/assembly"

# Reference
ref="${BASE_DIR}/ref/tifrunner_v2_filt.fa"

# Input assembly 
query="${BASE_DIR}/filt_assembly/TifNV-HG.asm.bp.p_ctg.25k.fa"
basename=$(basename "$query")
sample=${basename%.asm.bp.p_ctg.25k.fa}

# Output directories
OUT_FA="${BASE_DIR}/scaffold_25k"
OUT_RDS="${BASE_DIR}/rds_file"
OUT_PDF="${BASE_DIR}/pdf"

mkdir -p "$OUT_FA" "$OUT_RDS" "$OUT_PDF"

echo "Running Pteranodon for sample: $sample"
echo "Query:  $query"

# Parameters
SegLen=5000
MinQueryLen=1
threads=24


OUT_NAME="${sample}_SegLen${SegLen}_minQ${MinQueryLen}"

# Run Pteranodon
/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonBase.sh \
    -ref "$ref" \
    -query "$query" \
    -o "${OUT_FA}/${OUT_NAME}" \
    -SegLen $SegLen \
    -MinQueryLen $MinQueryLen \
    -auto 1 \
    -t $threads

# Move auxiliary outputs 
mv "${OUT_FA}/${OUT_NAME}"/*.rds "$OUT_RDS"/
mv "${OUT_FA}/${OUT_NAME}"/*.pdf "$OUT_PDF"/
