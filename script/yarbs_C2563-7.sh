#!/bin/bash
#SBATCH -J yarbs
#SBATCH -p khufu
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH -o /cluster/lab/clevenger/Ashmita/assembly/log/yarbs_%j.out
#SBATCH -e /cluster/lab/clevenger/Ashmita/assembly/log/yarbs_%j.err

module load cluster/minimap2/2.26
module load python/3.12.1-gcc-13.1.0

cd /cluster/home/aupadhyay/YARBS/python_scripts || exit 1

mkdir -p /cluster/lab/clevenger/Ashmita/assembly/yarbs

python minimap_prep.py \
  -r /cluster/lab/clevenger/Ashmita/assembly/ref/tifrunner_v2_filt.fa \
  -q /cluster/lab/clevenger/Ashmita/assembly/unfilt_assembly/C2563-7.asm.bp.p_ctg.fa \
  -o /cluster/lab/clevenger/Ashmita/assembly/yarbs/C2563-7_vs_tifrunner \
  -t 32 \
  --preset asm20 \
  --unique-length 10000 \
  --primary-only
