#!/bin/bash
#SBATCH -J ptero_TifHNV_autofix
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

# ---- INPUT ----
INPUT_DIR="/cluster/lab/clevenger/Ashmita/assembly/unfilt_assembly"
ref="/cluster/lab/clevenger/Ashmita/assembly/ref/tifrunner_v2_filt.fa"
query="${INPUT_DIR}/TifNV-HG.asm.bp.p_ctg.fa"

# ---- OUTPUT (SEPARATE DIR) ----
OUTDIR="/cluster/lab/clevenger/Ashmita/assembly/scaffold_autofix/TifHNV_manual_autofix"
mkdir -p "$OUTDIR"

# ---- HARD-CODED PREFIX ----
prefix="${OUTDIR}/TifHNV_Ptautofix_min01"

echo "Running Pteranodon"
echo "Query : $query"
echo "Ref   : $ref"
echo "Prefix: $prefix"

# Parameters
SegLen=1000
MinQueryLen=0.01
threads=32

# Run Pteranodon
/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonBase.sh \
    -ref "$ref" \
    -query "$query" \
    -o "$prefix" \
    -SegLen $SegLen \
    -MinQueryLen $MinQueryLen \
    -auto 0 \
    -t $threads

# Reformat FASTA to single-line sequences
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
  "${prefix}.fa" \
> "${prefix}.joined.fa"
