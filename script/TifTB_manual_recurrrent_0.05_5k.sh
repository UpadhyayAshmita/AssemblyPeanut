#!/bin/bash
#SBATCH -J ptero_TifTB_recurrent
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=60G
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_recurrent_%x_%j.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_recurrent_%x_%j.err"


# Load Khufu + Pteranodon env
source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnv.sh
source /cluster/projects/khufu/korani_projects/load_modules.sh

# --- paths from your directory listing ---
BASEDIR="/cluster/lab/clevenger/Ashmita/assembly/scaffold_autofix/TifTB_manual_autofix_min0.05_seg5k"
RUNDIR="${BASEDIR}/TifTB_Ptautofix_min01"

SCRIPT="${BASEDIR}/data-2026-02-02.sh"
QUERY="${RUNDIR}/Query.fa"
OUTFA="${BASEDIR}/TifTB_recurrent_out.fa"

cd "$RUNDIR"

/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonRecurrent.sh \
  -script "$SCRIPT" \
  -query "$QUERY" \
  -o "$OUTFA"

echo "DONE: $OUTFA"
