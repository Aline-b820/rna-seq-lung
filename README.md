# RNA-seq analysis of pulmonary response to *Toxoplasma gondii* infection

This repository contains an RNA-seq analysis of lung samples from wild-type (WT) and interferon receptor double knockout (DKO; Ifnar−/− × Ifngr−/−) mice infected with *Toxoplasma gondii*.  
The project was conducted as part of the Bioinformatics Master RNA-seq course and aims to reproduce and interpret key transcriptional signatures described in Singhania et al. (2019).

The analysis focuses on the effect of interferon receptor deficiency on the pulmonary immune response during infection.

---

## Experimental design

Lung samples were obtained from four experimental groups:
- WT Control  
- WT Case (infected with *T. gondii*)  
- DKO Control  
- DKO Case (infected with *T. gondii*)

The primary contrast of interest in this project is:
**Case DKO vs Case WT**, representing the effect of interferon receptor knockout during infection.

---

## Repository structure
```
rna-seq-lung/
├── data/
│   ├── fastq/              # Raw paired-end FASTQ files
│   └── metadata_lung.csv   # Sample metadata
│
├── refs/
│   ├── Mus_musculus.GRCm39.dna.primary_assembly.fa
│   ├── Mus_musculus.GRCm39.113.gtf
│   └── mm39_hisat2.*.ht2   # HISAT2 genome index files
│
├── qc/
│   ├── fastqc/             # FastQC reports for individual samples
│   └── multiqc/            # MultiQC summary report
│
├── align/
│   ├── *.bam               # Sorted BAM alignment files
│   ├── *.bam.bai           # BAM index files
│   ├── *_hisat2.log        # HISAT2 alignment logs
│   └── hisat2_summary_table.csv
│
├── counts/
│   ├── gene_counts.txt             # Gene-level counts (featureCounts)
│   └── gene_counts.txt.summary     # featureCounts assignment summary
│
├── scripts/
│   ├── 00_setup_project.sh
│   ├── 01_qc_fastqc.slurm
│   ├── 02_qc_multiqc.slurm
│   ├── 03_refs_mm39.slurm
│   ├── 04_align_mm39_per_sample.slurm
│   ├── 05_featureCounts_mm39.slurm
│   ├── locall_DESeq2_plots.R
│   ├── normalizad.R
│   ├── Expression-plots-for-Cxcl10-Ifng-Nos2.R
│   └── 07_GO_enrich_Case_DKO_vs_WT.R
│   
│   
│  
│
├── results/
│   ├── figures/
│   │   ├── PCA_all_samples_local.pdf
│   │   ├── Volcano_Case_DKO_vs_WT_local.pdf
│   │   └── GO_BP_Case_DKO_vs_WT_local.pdf
│   │
│   └── tables/
│       ├── DESeq2_Case_DKO_vs_WT_local.csv
│       ├── GO_BP_Case_DKO_vs_WT_local.csv
│       └──normalized_counts_Cxcl10_Ifng_Nos2_local.csv
│
│
│
├── README.md
└── LICENSE

```
## Analysis workflow

### 1. Quality control  
Raw paired-end RNA-seq reads were assessed using **FastQC (v0.11.9)** and summarized with **MultiQC (v1.11)** to evaluate read quality, adapter content, and sequencing depth.

### 2. Read alignment  
Reads were aligned to the *Mus musculus* reference genome (GRCm39) using **HISAT2 (v2.2.1)** with strand-specific RF settings.  
Alignments were converted to BAM format, sorted, and indexed using **SAMtools (v1.13)**.

### 3. Gene-level quantification  
Gene-level read counts were generated using **featureCounts** non–strand-specific (-s 0), based on Ensembl gene annotations.

### 4. Exploratory data analysis  
Differential expression analysis was performed in **R (v4.4.1)** using **DESeq2 (v1.46.0)**.  
Variance stabilizing transformation (VST) was applied for exploratory analyses, including principal component analysis (PCA), to assess sample clustering and data quality.

### 5. Differential expression analysis  
Differential expression was assessed using DESeq2.  
Genes with adjusted p-value (padj) < 0.05 and |log2 fold change| > 1 were considered significantly differentially expressed.

### 6. Functional enrichment analysis  
Gene Ontology (GO) over-representation analysis (Biological Process) was performed using **clusterProfiler (v4.14.6)** and the **org.Mm.eg.db** annotation package.

---

## Main results

Final results used in the report are located in:

### Figures
- `PCA_all_samples_local.pdf` – PCA of VST-transformed counts  
- `Volcano_Case_DKO_vs_WT_local.pdf` – Volcano plot for Case DKO vs WT  
- `GO_BP_Case_DKO_vs_WT_local.pdf` – GO Biological Process enrichment

### Tables
- `DESeq2_Case_DKO_vs_WT_local.csv` – Differential expression results  
- `normalized_counts_Cxcl10_Ifng_Nos2_local.csv` – Normalized counts for selected interferon-related genes  
- `GO_BP_Case_DKO_vs_WT_local.csv` – GO enrichment results

---

## Reproducibility

All analyses were performed using documented scripts and fixed software versions.  
The repository is structured to allow another bioinformatician to reproduce the analysis starting from raw FASTQ files through to the final figures and tables.

---

## References

Singhania A, Graham CM, Gabrysova L, et al.  
Transcriptional profiling unveils type I and II interferon networks in *Toxoplasma gondii* infection.  
*Nature Communications*. 2019;10:2887.

Love MI, Huber W, Anders S.  
Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.  
*Genome Biology*. 2014;15:550.

Kim D, Paggi JM, Park C, Bennett C, Salzberg SL.  
Graph-based genome alignment and genotyping with HISAT2.  
*Nature Biotechnology*. 2019;37:907–915.

Liao Y, Smyth GK, Shi W.  
featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features.  
*Bioinformatics*. 2014;30:923–930.

Wu T, Hu E, Xu S, et al.  
clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.  
*The Innovation*. 2021;2:100141.
