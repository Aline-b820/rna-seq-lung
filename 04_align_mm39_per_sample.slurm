#!/usr/bin/env bash
#SBATCH -J align_mm39
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH -o align/align_mm39_%A_%a.out
#SBATCH -e align/align_mm39_%A_%a.err

set -euo pipefail

module load HISAT2
module load SAMtools

PROJ=/data/users/abarbosa/projects/rna-seq-lung
FASTQ=$PROJ/data/fastq
ALIGN=$PROJ/align
REFS=$PROJ/refs
IDX=$REFS/mm39_hisat2

mkdir -p "$ALIGN"

# lista de samples (em ORDEM!)
SAMPLES=(
  SRR7821918
  SRR7821919
  SRR7821920
  SRR7821921
  SRR7821922
  SRR7821923
  SRR7821924
  SRR7821925
  SRR7821927
  SRR7821937
  SRR7821938
  SRR7821939
  SRR7821940
  SRR7821941
  SRR7821942
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1=$FASTQ/${SAMPLE}_1.fastq.gz
R2=$FASTQ/${SAMPLE}_2.fastq.gz

echo ">> alinhando $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"
echo "Índice: $IDX"

# If you already have a BAM from the previous round, don't overwrite it accidentally.
if [[ -s "$ALIGN/${SAMPLE}.bam" ]]; then
  echo ">> WARNING: $ALIGN/${SAMPLE}.bam already exists. Delete or rename it before running again."
  exit 1
fi

time hisat2 -x "$IDX" \
  -1 "$R1" -2 "$R2" \
  --dta -p 8 \
  --summary-file "$ALIGN/${SAMPLE}_hisat2_summary.txt" \
  2> "$ALIGN/${SAMPLE}_hisat2.log" \
| samtools sort -@ 8 -o "$ALIGN/${SAMPLE}.bam"

samtools index "$ALIGN/${SAMPLE}.bam"

echo ">> alinhamento concluído para $SAMPLE"

How to run:

#sbatch --array=0-14 scripts/02_align_mm39_per_sample.slurm
#Note: Adjust the array according to the number of samples in the SAMPLES list above; it will create a job for each sample.
#Each job will take the sample corresponding to the array index, each with a different SAMPLE,
#and generate the corresponding batch file; That is, if there are 15 samples, the array will go from 0-14 (15 jobs).
#And each job will generate a different batch file, one per sample.

#####

#Check progress with:
#squeue -u $USER -n align_mm39
#####
#View log of a sample (e.g., the first in the array, TASK_ID=0):
#ls -lh align/align_mm39_*_0.out
#cat align/align_mm39_*_0.out
#cat align/align_mm39_*_0.err
