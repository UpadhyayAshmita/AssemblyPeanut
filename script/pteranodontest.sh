#!/bin/bash
#SBATCH -J ptero_AU23_10
#SBATCH --time=24:00:00
#SBATCH -c 24
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=180G
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%J.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/ptero_%x_%J.err"


#   Load Khufu + Pteranodon env
source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnv.sh
source "$khufu_dir"/utilities/load_modules.sh

#   Input paths
ref="/cluster/lab/clevenger/Ashmita/assembly/ref/tifrunner_v2.fa"
query="/cluster/lab/clevenger/Ashmita/assembly/AU-23-10.asm.bp.p_ctg.50k.fa"
out="/cluster/lab/clevenger/Ashmita/assembly/AU-23-10_scaffold"

# Parameters
SegLen=500000       # segment length to chop reference
MinQueryLen=1000    # minimum contig length
threads=24
#   Run Pteranodon
/cluster/projects/khufu/korani_projects/Pteranodon/scripts/PteranodonBase.sh \
    -ref $ref \
    -query $query \
    -o $out \
    -SegLen $SegLen \
    -MinQueryLen $MinQueryLen \
    -t $threads
