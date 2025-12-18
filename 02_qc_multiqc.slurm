#!/usr/bin/env bash
#SBATCH -J qc_multiqc_lung
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH -o qc/qc_multiqc_%j.out
#SBATCH -e qc/qc_multiqc_%j.err

set -euo pipefail

PROJ=/data/users/abarbosa/projects/rna-seq-lung
QC_DIR=$PROJ/qc
FQ_OUT=$QC_DIR/fastqc
MQ_OUT=$QC_DIR/multiqc

mkdir -p "$MQ_OUT"


module load MultiQC 2>/dev/null || echo "Warning: Adjust to use apptainer if there is no MultiQC module."

echo ">> Rodando MultiQC"
echo ">> Entrada: $FQ_OUT, $PROJ/align, $PROJ/counts"
echo ">> Saída:   $MQ_OUT"

multiqc \
  "$FQ_OUT" \
  "$PROJ/align" \
  "$PROJ/counts" \
  -o "$MQ_OUT"

echo ">> MultiQC completed."

