library(tidyverse)

# Load normalized counts
df <- read.csv("normalized_counts_all_samples.csv")

# Order groups for consistency
df$group <- factor(df$group, levels = c(
  "Lung_WT_Control", "Lung_WT_Case",
  "Lung_DKO_Control", "Lung_DKO_Case"))

genes_of_interest <- c("Cxcl10", "Ifng", "Nos2")

df_long <- df %>%
  filter(gene %in% genes_of_interest) %>%
  pivot_longer(cols = starts_with("Lung_"), 
               names_to = "group", values_to = "count")

# Clean group names
df_long$group <- df_long$group %>% 
  str_replace("Lung_", "") %>%
  factor(levels = c("WT_Control", "WT_Case", "DKO_Control", "DKO_Case"))

# Plot function
plot_gene <- function(g) {
  ggplot(df_long %>% filter(gene == g),
         aes(x = group, y = count, fill = group)) +
    geom_bar(stat = "summary", fun = "mean") +
    geom_point(size = 3, alpha = 0.6) +
    theme_bw(base_size = 15) +
    ggtitle(paste("Expression of", g)) +
    ylab("Normalized DESeq2 count") +
    xlab("") +
    theme(legend.position = "none")
}

# Save plots
for (g in genes_of_interest) {
  ggsave(paste0("Expression_", g, ".pdf"),
         plot = plot_gene(g), width = 6, height = 5)
}
