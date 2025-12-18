#!/usr/bin/env bash

echo "=========================================="
echo "RNA-seq project setup"
echo "Creating project directories..."
echo "=========================================="

# PROJ = directory where the script is being executed
PROJ=/data/users/abarbosa/projects/rna-seq-lung

mkdir -p $PROJ/data/fastq
mkdir -p $PROJ/data
mkdir -p $PROJ/refs
mkdir -p $PROJ/align
mkdir -p $PROJ/counts
mkdir -p $PROJ/results/figures
mkdir -p $PROJ/results/tables
mkdir -p $PROJ/scripts
mkdir -p $PROJ/qc

echo "Structure created:"
tree -L 3 $PROJ

echo "=========================================="
echo "Setup complete!"
echo "=========================================="

