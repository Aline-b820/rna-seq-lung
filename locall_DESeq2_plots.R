############################################################
## local_DESeq2_plots.R
## Local analysis in RStudio based on:
##   - gene_counts.txt
##   - metadata_lung.csv
############################################################

## 0)Pactes --------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

theme_set(theme_bw(base_size = 13))

## 1) Pathsand data reading ---------------------------------------------

proj_dir    <- getwd()
counts_file <- file.path(proj_dir, "gene_counts.txt")
meta_file   <- file.path(proj_dir, "metadata_lung.csv")

cat(">> Using project in:", proj_dir, "\n")
cat(">> Counts:", counts_file, "\n")
cat(">> Metadata:", meta_file, "\n\n")

## 1.1) Read counts from featureCounts ------------------------------------
cts <- read.delim(counts_file,
                  comment.char = "#",
                  check.names  = FALSE)

## First column = Geneid; 2–6 are annotations
rownames(cts) <- cts$Geneid
cts <- cts[, -(1:6)]

##Column names → SRRxxxx
samples <- sub("\\.bam$", "", basename(colnames(cts)))
colnames(cts) <- samples

cat(">> Samples in the counting columns:\n")
print(samples)

## 1.2) read metadata -------------------------------------------------------
meta <- read_csv(meta_file, show_col_types = FALSE)
meta <- as.data.frame(meta)

## Sort metadata in the same order as the count columns.
meta <- meta[match(samples, meta$sample), ]
rownames(meta) <- meta$sample

cat("\n>> Metadata aligned with counts:\n")
print(meta[, c("sample","genotype","condition","group")])

meta$genotype  <- factor(meta$genotype,  levels = c("WT","DKO"))
meta$condition <- factor(meta$condition, levels = c("Control","Case"))
meta$group     <- factor(meta$group)

## 2)To create DESeqDataSet + DESeq + VST ------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(cts),
  colData   = meta,
  design    = ~ group
)

## Low count filter
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

cat("\n>> Number of genes after filtering (soma counts >= 10):", nrow(dds), "\n\n")

dds <- DESeq(dds)

## VST to PCA (blind=TRUE)
vsd     <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

## 3) PCA exploratory --------------------------------------------------

cat(">> Generating PCA (vst) with plotPCA... \n")

pcaData <- plotPCA(vsd,
                   intgroup = c("genotype","condition"),
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData,
                aes(x = PC1, y = PC2,
                    color = genotype,
                    shape = condition,
                    label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("WT" = "steelblue4", "DKO" = "firebrick3")) +
  scale_shape_manual(values = c("Control" = 16, "Case" = 17)) +
  labs(
    title = "PCA of VST-transformed counts (all samples)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Genotype",
    shape = "Condition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

dir.create("results/figures_local", showWarnings = FALSE, recursive = TRUE)
pca_pdf <- file.path("results/figures_local", "PCA_all_samples_local.pdf")
pca_png <- file.path("results/figures_local", "PCA_all_samples_local.png")

ggsave(pca_pdf, p_pca, width = 6, height = 5)
ggsave(pca_png, p_pca, width = 6, height = 5, dpi = 300)

cat("  -> PCA salvo em:\n     -", pca_pdf, "\n     -", pca_png, "\n\n")

## 4) DE for main contrast: Case DKO vs WT -------------------------
cat(">> DE for main contrast: Case DKO vs WT (Lung_DKO_Case vs Lung_WT_Case)...\n")
res_case <- results(dds,
                    contrast = c("group", "Lung_DKO_Case", "Lung_WT_Case"))

res_case_df <- as.data.frame(res_case)
res_case_df$gene_id <- rownames(res_case_df)

## liguar gene_id, gene symbol via GTF annotation;
## Ensembl ID.
dir.create("results/tables_local", showWarnings = FALSE, recursive = TRUE)
write.csv(res_case_df,
          "results/tables_local/DESeq2_Case_DKO_vs_WT_local.csv",
          row.names = FALSE)

padj_cut <- 0.05
lfc_cut  <- 1

sig_case <- res_case_df %>%
  filter(!is.na(padj) & padj < padj_cut)

case_up   <- sum(sig_case$log2FoldChange >  lfc_cut, na.rm = TRUE)
case_down <- sum(sig_case$log2FoldChange < -lfc_cut, na.rm = TRUE)

cat("   DEGs (padj <", padj_cut, ", |log2FC| >", lfc_cut, "):\n")
cat("      UP   =", case_up, "\n")
cat("      DOWN =", case_down, "\n\n")

## 5) Volcano main contrast plot ----------
cat(">> Gerando volcano (Case DKO vs WT)...\n")

plot_volcano <- function(df,
                         title,
                         outfile,
                         lfc_cut = 1,
                         padj_cut = 0.05,
                         highlight_ids = c("ENSMUSG00000034855",  # Cxcl10
                                           "ENSMUSG00000055170",  # Ifng
                                           "ENSMUSG00000020826")) { # Nos2
  
  df2 <- df %>%
    filter(!is.na(padj)) %>%
    mutate(
      logFC    = log2FoldChange,
      neglog10 = -log10(padj),
      signif   = padj < padj_cut & abs(logFC) > lfc_cut
    )
  
  ## labels: genes of interest + top 5 up/down by padj  top_up <- df2 %>%
    filter(signif, logFC > 0) %>%
    arrange(padj) %>%
    head(5) %>%
    pull(gene_id)
  
  top_down <- df2 %>%
    filter(signif, logFC < 0) %>%
    arrange(padj) %>%
    head(5) %>%
    pull(gene_id)
  
  label_ids <- unique(c(highlight_ids, top_up, top_down))
  df2$label <- ifelse(df2$gene_id %in% label_ids, df2$gene_id, NA)
  
  p <- ggplot(df2, aes(x = logFC, y = neglog10)) +
    geom_point(aes(color = signif), alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"),
                       guide = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_cut),
               linetype = "dashed", color = "grey50") +
    geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.3,
      point.padding = 0.1,
      min.segment.length = 0
    ) +
    labs(
      title = title,
      x = "log2 fold change (DKO vs WT, Case)",
      y = "-log10(FDR)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
  
  ggsave(outfile, p, width = 6, height = 5)
  cat("   -> Volcano salvo em:", outfile, "\n")
}

volc_pdf <- file.path("results/figures_local", "Volcano_Case_DKO_vs_WT_local.pdf")
plot_volcano(res_case_df,
             title   = "Volcano – Case DKO vs WT",
             outfile = volc_pdf,
             lfc_cut = lfc_cut,
             padj_cut = padj_cut)

cat("\n")

## 6) Normalized counts of Cxcl10, Ifng, Nos2

cat(">> Extracting normalized counts from Cxcl10, Ifng, Nos2...\n")
norm_counts <- counts(dds, normalized = TRUE)

ids_interest <- c("ENSMUSG00000034855",  # Cxcl10
                  "ENSMUSG00000055170",  # Ifng
                  "ENSMUSG00000020826")  # Nos2
symbols_interest <- c("Cxcl10","Ifng","Nos2")
names(symbols_interest) <- ids_interest

sub_counts <- norm_counts[ids_interest, , drop = FALSE]

norm_tab <- as.data.frame(t(sub_counts))
norm_tab$sample <- rownames(norm_tab)

norm_tab <- norm_tab %>%
  left_join(meta, by = c("sample" = "sample")) %>%
  relocate(sample, genotype, condition, group)

## Rename gene columns
colnames(norm_tab)[colnames(norm_tab) %in% ids_interest] <- symbols_interest

out_nc <- "results/tables_local/normalized_counts_Cxcl10_Ifng_Nos2_local.csv"
write.csv(norm_tab, out_nc, row.names = FALSE)

cat(" -> Normalized count table saved in:\n ", out_nc, "\n")

cat("\n=== END OF LOCAL SCRIPT (PCA + Volcano + 3 genes) ===\n")
