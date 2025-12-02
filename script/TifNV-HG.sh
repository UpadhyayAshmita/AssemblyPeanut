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

INPUT_DIR="/cluster/lab/clevenger/Ashmita/assembly"
ref="${INPUT_DIR}/ref/tifrunner_v2.fa"

query="${INPUT_DIR}/TifNV-HG.asm.bp.p_ctg.50k.fa"
basename=$(basename "$query")
sample=${basename%.asm.bp.p_ctg.50k.fa}

# Output directory
out="${INPUT_DIR}/${sample}_scaffold"
mkdir -p "$out"

echo "Running Pteranodon for sample: $sample"
echo "Query:  $query"
echo "Output: $out"

# Parameters
SegLen=1000
MinQueryLen=1
threads=24

# Run Pteranodon
/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonBase.sh \
    -ref "$ref" \
    -query "$query" \
    -o "$out" \
    -SegLen $SegLen \
    -MinQueryLen $MinQueryLen \
    -auto 1 \
    -t $threads
