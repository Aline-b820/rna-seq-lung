## 13_normalized_counts.R
## Generates normalized counts and extracts Cxcl10, Ifng, Nos2
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

proj_dir    <- getwd()
counts_file <- file.path(proj_dir, "gene_counts.txt")
meta_file   <- file.path(proj_dir, "metadata_lung.csv")

cat(">> Reading counts of:\n   ", counts_file, "\n")
cts <- read.delim(counts_file,
                  comment.char = "#",
                  check.names  = FALSE)

## The first column is Geneid, then Chr, Start, End, Strand, Length.
rownames(cts) <- cts$Geneid
cts <- cts[, -(1:6)]

## colnames: full paths of the BAMs in the cluster -> becomes SRRxxxx
samples <- sub("\\.bam$", "", basename(colnames(cts)))
colnames(cts) <- samples

cat("Samples detectados nas contagens:\n")
print(samples)

cat(">> Read metadata de:\n   ", meta_file, "\n")
meta <- read_csv(meta_file, show_col_types = FALSE) |>
  as.data.frame()

## align metadata with count columns
meta <- meta[match(samples, meta$sample), ]
rownames(meta) <- meta$sample

if (any(is.na(rownames(meta)))) {
  stop("Some sample counts were not found in the metadata!")
}

cat("Metadata aligned with counts:\n")
print(meta[, c("sample", "genotype", "condition", "group")])

## factors with ordered levels
meta$genotype  <- factor(meta$genotype,  levels = c("WT", "DKO"))
meta$condition <- factor(meta$condition, levels = c("Control", "Case"))

## ----------------------------------------------------------
##DESeq2: for normalization (size factors) only.
## ----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(cts),
  colData   = meta,
  design    = ~ genotype + condition + genotype:condition
)

## filter out very low-level genes (as we did before)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

out_norm <- file.path(proj_dir, "normalized_counts_all_samples.csv")
write.csv(as.data.frame(norm_counts),
          out_norm,
          quote = FALSE)

cat(">> Normalized counts saved in:\n ", out_norm, "\n")

## ----------------------------------------------------------
## 3 genes of interest: Cxcl10, Ifng, Nos2
## (Mus musculus Ensembl IDs)
## ----------------------------------------------------------
genes_interest <- c(
  Cxcl10 = "ENSMUSG00000034855",
  Ifng   = "ENSMUSG00000055170",
  Nos2   = "ENSMUSG00000020826"
)

present <- genes_interest[genes_interest %in% rownames(norm_counts)]
if (length(present) == 0) {
  stop("None of the genes of interest were found in gene_counts.txt!")
}

go_tab <- norm_counts[unname(present), , drop = FALSE]
go_tab <- as.data.frame(go_tab)
go_tab$gene_id <- rownames(go_tab)
go_tab$gene    <- names(present)[match(go_tab$gene_id, present)]
go_tab <- go_tab[, c("gene", "gene_id", samples)]

##Calculate average per group (Lung_WT_Case, etc.)
meta$group <- factor(meta$group)

group_means <- sapply(levels(meta$group), function(g) {
  idx <- meta$group == g
  rowMeans(go_tab[, samples[idx], drop = FALSE])
})

group_means <- as.data.frame(group_means)
group_means$gene    <- go_tab$gene
group_means$gene_id <- go_tab$gene_id
group_means <- group_means[, c("gene", "gene_id", levels(meta$group))]

out_go <- file.path(proj_dir, "genes_of_interest_normalized_counts.csv")
write.csv(group_means, out_go, row.names = FALSE)

cat(">> Tabela de Cxcl10 / Ifng / Nos2 por grupo salva em:\n   ",
    out_go, "\n")

cat("=== Script 13_normalized_counts.R finished. ===\n")
