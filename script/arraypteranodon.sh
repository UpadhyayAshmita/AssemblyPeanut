#!/bin/bash
#SBATCH -J ptero_all_fa
#SBATCH --time=48:00:00
#SBATCH -c 24
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=180G
#SBATCH --array=1-12          
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%A_%a.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%A_%a.err"

# Load Khufu + Pteranodon env
source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnv.sh
source /cluster/projects/khufu/korani_projects/load_modules.sh

# Input directory

INPUT_DIR="/cluster/lab/clevenger/Ashmita/assembly"
ref="${INPUT_DIR}/ref/tifrunner_v2.fa"

files=(${INPUT_DIR}/*50k.fa)

query="${files[$SLURM_ARRAY_TASK_ID-1]}"
basename=$(basename "$query")
sample=${basename%.asm.bp.p_ctg.50k.fa}

# Output directory
out="${INPUT_DIR}/${sample}_scaffold"
mkdir -p "$out"

echo " Running Pteranodon for sample: $sample"
echo " Query: $query"
echo " Output: $out"


# Pteranodon parameters
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
