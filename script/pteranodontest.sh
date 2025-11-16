#!/bin/bash
#SBATCH -J pteranodon_test
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=128G
#SBATCH -o "/cluster/lab/clevenger/Ashmita/assembly/log/pteranodon_%x_%J.out"
#SBATCH -e "/cluster/lab/clevenger/Ashmita/assembly/log/pteranodon_%x_%J.err"

# Path to Pteranodon installation
pteranodon="/cluster/home/wkorani/korani_aps/pteranodon"
t=8

# Reference genome (Tifrunner V2)
ref="/cluster/lab/clevenger/Ashmita/assembly/ref/tifrunner_v2.fa"

# One of the 50k filtered contigs
query="/cluster/lab/clevenger/Ashmita/assembly/AU-23-10.asm.bp.p_ctg.50k.fa"

# Output folder
out="/cluster/lab/clevenger/Ashmita/assembly/AU-23-10_scaffold"

# Minimum contig length threshold
len=1000

# Minimum percent identity
min=1

# Run Pteranodon
"$pteranodon"/get_plot_based_on_local_alignment.sh "$pteranodon" "$t" "$ref" "$query" "$out" "$len" "$min"
