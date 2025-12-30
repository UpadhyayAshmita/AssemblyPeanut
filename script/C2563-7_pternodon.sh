#!/bin/bash
#SBATCH -J ptero_C2563-7_autofix
#SBATCH --time=48:00:00
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=180G
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_autofix_%x_%j.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_autofix_%x_%j.err"

# Load Khufu + Pteranodon env
source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnv.sh
source /cluster/projects/khufu/korani_projects/load_modules.sh

INPUT_DIR="/cluster/lab/clevenger/Ashmita/assembly/unfilt_assembly"
ref="${INPUT_DIR}/ref/tifrunner_v2_filt.fa"

query="${INPUT_DIR}/C2563-7.asm.bp.p_ctg.fa"
basename=$(basename "$query")
sample=${basename%.asm.bp.p_ctg.fa}

# Output directory
out="${INPUT_DIR}/${sample}_scaffold_autofix_min.01"
mkdir -p "$out"

echo "Running Pteranodon for sample: $sample"
echo "Query:  $query"
echo "Output: $out"

# Parameters
SegLen=1000
MinQueryLen=.01
threads=32

# Run Pteranodon
/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonBase.sh \
    -ref "$ref" \
    -query "$query" \
    -o "$out" \
    -SegLen $SegLen \
    -MinQueryLen $MinQueryLen \
    -auto 1 \
    -t $threads



awk '/^>/ {
        if (seq) print seq;
        print;
        seq="";
        next
     }
     {
        seq = seq $0
     }
     END {
        print seq
     }' \
  "${out}.fa" \
> "${out}.joined.fa"
