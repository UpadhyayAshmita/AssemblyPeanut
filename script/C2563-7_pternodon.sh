#!/bin/bash
#SBATCH -J ptero_C2563-7
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
ref="${INPUT_DIR}/ref/tifrunner_v2_filt.fa"

query="${INPUT_DIR}/C2563-7.asm.bp.p_ctg.20k.fa"
basename=$(basename "$query")
sample=${basename%.asm.bp.p_ctg.20k.fa}

# Output directory
out="${INPUT_DIR}/${sample}_scaffold_20k_min0.05"
mkdir -p "$out"

echo "Running Pteranodon for sample: $sample"
echo "Query:  $query"
echo "Output: $out"

# Parameters
SegLen=5000
MinQueryLen=0.05
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
