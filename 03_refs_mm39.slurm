#!/usr/bin/env bash
#SBATCH -J idx_mm39
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH -o refs/idx_mm39_%j.out
#SBATCH -e refs/idx_mm39_%j.err

set -euo pipefail

PROJ=/data/users/abarbosa/projects/rna-seq-lung
REFS=$PROJ/refs
# load HISAT2 from cluster module
module load HISAT2

FA=$REFS/Mus_musculus.GRCm39.dna.primary_assembly.fa

echo ">> usando FASTA: $FA"
head -n 3 "$FA"

echo ">> construindo índice HISAT2 mm39_hisat2 ..."
hisat2-build -p 8 "$FA" "$REFS/mm39_hisat2"

echo ">> índice criado em $REFS/mm39_hisat2.*.ht2"
