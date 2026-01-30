#!/usr/bin/env bash
#SBATCH -J qc_fastqc_lung
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH -o qc/qc_fastqc_%j.out
#SBATCH -e qc/qc_fastqc_%j.err

set -euo pipefail

# Project Path
PROJ=/data/users/abarbosa/projects/rna-seq-lung

# Directories
FASTQ=$PROJ/data/fastq
QC_DIR=$PROJ/qc
FQ_OUT=$QC_DIR/fastqc

# Container FastQC path
FASTQC_SIF=/containers/apptainer/fastqc-0.12.1.sif

mkdir -p "$FQ_OUT"

echo ">> Rodando FastQC em todos os FASTQ de $FASTQ"
echo ">> Resultados serão salvos em $FQ_OUT"

for fq in "$FASTQ"/*.fastq.gz; do
  echo ">> FastQC: $fq"
  apptainer exec "$FASTQC_SIF" \
    fastqc "$fq" -o "$FQ_OUT"
done

echo ">> FastQC concluído."

