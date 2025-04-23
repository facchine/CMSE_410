library(DESeq2)
library(IRanges)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(textshaping)
library(Rcpp)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggrepel)

counts_matrix <- read.table("path/to/star_salmon/salmon.merged.gene_counts.tsv", header=TRUE, row.names=1, sep="\t")
print(dim(counts_matrix))   # Check dimensions (rows=genes, cols=samples+metadata?)

g2s <- data.frame(
  gene_id = rownames(counts_matrix), 
  gene_symbol = counts_matrix[, 1]
  )

counts_matrix <- counts_matrix[, -1]

counts_matrix <- as.matrix(counts_matrix)
View(counts_matrix)

counts_matrix_rounded <- round(counts_matrix)
View(counts_matrix_rounded)

counts_filtered <- counts_matrix_rounded[rowSums(counts_matrix_rounded) > 1, ]

save(counts_filtered, counts_matrix, g2s, file = "results/count_files.RData", verbose = T)

counts_matrix <- counts_matrix %>%
  select(-"gene_name")
  
deseq_samples <- data.frame(
  sample_id = colnames(counts_matrix)
)

split_values <- strsplit(deseq_samples$sample_id, "_")

samplesheet <- read.csv("/path/to/deseqsamplesheet.csv")
#samplesheet <- samplesheet[, -1]

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = samplesheet,
                              design = ~ time_point)

dds <- DESeq(dds)

resultsNames(dds)

rlog_counts <- rlog(dds, blind = TRUE)

rlog_counts_matrix <- assay(rlog_counts)

write_rds(rlog_counts_matrix, "results/rlog_counts_yeast.rds")

res1 <- results(dds, contrast = c("time_point","log", "1"))

res1_df <- data.frame(time = "log", res1)
view(res1_df)

res1_df <- rownames_to_column(res1_df, "gene_id")
View(res1_df)

sum(is.na(res1_df$padj)) / nrow(res1_df)

filtered_res1_df <- res1_df %>%
  filter(padj < 0.05)

nrow(filtered_res1_df)

filtered_res1_df_2 <- res1_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

nrow(filtered_res1_df_2)

genes_of_interest <- cat(paste(filtered_res1_df_2$gene_id, collapse = "\n"))

genes_of_interest_log <- as.data.frame(filtered_res1_df_2$gene_id, collapse = "\n")


g2s <- tibble(gene_id = rownames(counts_matrix), gene_symbol = counts_matrix[, 1])
counts_matrix <- counts_matrix[, -1]

samplesheet$strain <- factor(samplesheet$strain)
samplesheet$time_point <- factor(samplesheet$time_point, levels = c("log", "1", "3", "5"))

valid_strains <- samplesheet %>%
  filter(time_point %in% c("log", "1")) %>%
  group_by(strain, time_point) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(strain) %>%
  summarise(timepoints_present = n(), .groups = "drop") %>%
  filter(timepoints_present == 2) %>%
  pull(strain)


dir.create("results/strainwise_comparisons", recursive = TRUE, showWarnings = FALSE)

for (strain in valid_strains) {
  message("Processing: ", strain)
  
  ss_sub <- samplesheet %>%
    filter(strain == !!strain, time_point %in% c("log", "1"))
  
  cm_sub <- counts_filtered[, ss_sub$sample_id]
  rownames(ss_sub) <- ss_sub$sample_id
  
  
  cm_sub <- counts_filtered[, rownames(ss_sub)]
  
  dds_sub <- DESeqDataSetFromMatrix(countData = cm_sub, colData = ss_sub, design = ~ time_point)
  dds_sub <- DESeq(dds_sub)
  
  res <- results(dds_sub, contrast = c("time_point", "log", "1"))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(g2s, by = "gene_id")
  
  write_csv(res_df, paste0("results/strainwise_comparisons/", strain, "_log_vs_day1_results.csv"))
  
  sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  message(strain, ": ", nrow(sig_genes), " DEGs with padj < 0.05 and |log2FC| > 1")
  
  volcano_data <- res_df %>%
    mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1)
  
  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(
      data = subset(volcano_data, sig),
      aes(label = gene_id),
      size = 3, max.overlaps = 20
    ) +
    theme_minimal() +
    labs(
      title = paste0("Volcano: ", strain, " (log vs day 1)"),
      x = "Log2 Fold Change",
      y = "-log10(padj)"
    )
  
  ggsave(
    filename = paste0("results/strainwise_comparisons/", strain, "_volcano.pdf"),
    plot = volcano_plot,
    width = 8, height = 6
  )
  
  pdf(paste0("results/strainwise_comparisons/", strain, "_MAplot.pdf"))
  plotMA(res, main = paste0("MA plot: ", strain), ylim = c(-5, 5))
  dev.off()
  
  top_genes <- sig_genes %>% top_n(30, wt = abs(log2FoldChange)) %>% pull(gene_id)
  rlog_mat <- assay(rlog(dds_sub, blind = TRUE))
  
  pheatmap(rlog_mat[top_genes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           annotation_col = ss_sub,
           main = paste0("Top 30 DEGs: ", strain),
           filename = paste0("results/strainwise_comparisons/", strain, "_heatmap.pdf"))
}

de_results <- read_csv("results/strainwise_comparisons/yMK8393_log_vs_day1_results.csv")
top_genes <- de_results %>% filter(padj < 0.05, abs(log2FoldChange) > 1)

de_results <- read_csv("results/strainwise_comparisons/yXL006_log_vs_day1_results.csv")
top_genes <- de_results %>% filter(padj < 0.05, abs(log2FoldChange) > 1)

de_results <- read_csv("results/strainwise_comparisons/yXL171_log_vs_day1_results.csv")
top_genes <- de_results %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
