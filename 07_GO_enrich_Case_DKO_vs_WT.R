#########################################################
## 07_GO_enrich_Case_DKO_vs_WT.R
## GO (BP) enrichment analysis for the primary contrast agent:
## Case DKO vs WT (lung)
############################################################

## 0) Load packages ----
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

####b1) Paths 
proj_dir <- "C:/Users/aline/RNA sec"

de_file  <- file.path(
  proj_dir,
  "results", "tables_local",
  "DESeq2_Case_DKO_vs_WT_local.csv"
)

out_go_csv <- file.path(
  proj_dir,
  "results", "tables_local",
  "GO_BP_Case_DKO_vs_WT_local.csv"
)

out_go_pdf <- file.path(
  proj_dir,
  "results", "figures_local",
  "GO_BP_Case_DKO_vs_WT_local.pdf"
)

cat(">> read DE de:\n   ", de_file, "\n\n")

res <- read.csv(de_file, stringsAsFactors = FALSE)

## 2) Definir corte de DE ----
padj_cut <- 0.05
lfc_cut  <- 1

sig_genes <- res %>%
  filter(!is.na(padj),
         padj < padj_cut,
         abs(log2FoldChange) > lfc_cut)

cat("Genes DE (padj <", padj_cut,
    "e |log2FC| >", lfc_cut, "):",
    nrow(sig_genes), "\n")

## 3) Enrichment Gene List (ENSEMBL) ----
gene_DE      <- unique(sig_genes$gene_id)   # DE genes
gene_universe <- unique(res$gene_id)        # todos os genes testados

length(gene_DE)
length(gene_universe)

## 4) enrichGO (BP) ----
ego_bp <- enrichGO(
  gene          = gene_DE,
  universe      = gene_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",          # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE           # converts to gene symbols
)

## 5) View result in the console ----
ego_bp
head(as.data.frame(ego_bp))[ , 1:6]

## 6) Save table ----
dir.create(file.path(proj_dir, "results", "tables_local"),
           showWarnings = FALSE, recursive = TRUE)

write.csv(
  as.data.frame(ego_bp),
  out_go_csv,
  row.names = FALSE
)

cat(">> Tabela GO (BP) salva em:\n   ", out_go_csv, "\n")

## 7) Save figure (dotplot) -
dir.create(file.path(proj_dir, "results", "figures_local"),
           showWarnings = FALSE, recursive = TRUE)

pdf(out_go_pdf, width = 7, height = 6)
print(dotplot(ego_bp, showCategory = 15,
              title = "GO BP enrichment – Case DKO vs WT (lung)"))
dev.off()

cat(">> Figure GO (BP) saved in:\n   ", out_go_pdf, "\n\n")

cat("=== GO enrichment (BP) completed successfully ===\n")
