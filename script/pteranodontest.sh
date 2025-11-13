#!/bin/bash
#SBATCH -J pteranodon_test
#SBATCH --time=12:00:00
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -p khufu
#SBATCH --mem=128G
#SBATCH -o "stds/stdout_%x_%J"
#SBATCH -e "stds/stderr_%x_%J"

# Path to Pteranodon installation
pteranodon="/cluster/home/wkorani/korani_aps/pteranodon"
t=8

# Reference genome (Tifrunner V2)
ref="/cluster/lab/clevenger/Ashmita/assembly/ref/tifrunner_v2.fa"

# One of your query contigs (50k filtered)
query="/cluster/lab/clevenger/Ashmita/assembly/TifTB.asm.bp.p_ctg.50k.fa"

# Output folder name
out="TifTB_TR2_50k_scaffold "

# Minimum contig length threshold (1000 bp)
len=1000

# Minimum percent identity
min=1

# Run Pteranodon

"$pteranodon"/get_plot_based_on_local_alignment.sh "$pteranodon" "$t" "$ref" "$query" "$out" "$len" "$min"
