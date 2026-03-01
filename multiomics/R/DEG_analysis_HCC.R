# DEG analysis - hepatocellular carcinoma dataset
# using TCGA data downloaded from GDC portal
# DM - started this around Jan 2023, kept updating

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

# load the count matrix
# NOTE: make sure row names are gene symbols not ensembl IDs
# i had to convert them separately using biomaRt - see convert_ids.R

counts <- read.csv("data/hcc_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("data/hcc_metadata.csv", row.names = 1)

# quick sanity check - always do this first
dim(counts)
all(colnames(counts) == rownames(meta))  # must be TRUE

# some samples had very low library sizes in this dataset
# removing them before analysis
lib_sizes <- colSums(counts)
hist(log10(lib_sizes), breaks = 30, main = "Library size distribution")

low_libs <- names(lib_sizes[lib_sizes < 500000])
if(length(low_libs) > 0){
  message("removing ", length(low_libs), " low quality samples")
  counts <- counts[, !colnames(counts) %in% low_libs]
  meta   <- meta[!rownames(meta) %in% low_libs, ]
}

# filter lowly expressed genes
# keeping genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, ]
message(nrow(counts) - nrow(counts_filt), " genes removed after filtering")

# build dds object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData   = meta,
  design    = ~ condition
)

# set reference level - control/normal goes first
dds$condition <- relevel(dds$condition, ref = "normal")

# run DESeq2
dds <- DESeq(dds)

# extract results - tumor vs normal
res <- results(dds,
               contrast      = c("condition", "tumor", "normal"),
               alpha         = 0.05,
               pAdjustMethod = "BH")

summary(res)

# lfc shrinkage - apeglm is better than normal shrinkage for most cases
# ref: Zhu et al 2018
resultsNames(dds)
res_shrunk <- lfcShrink(dds,
                         coef = "condition_tumor_vs_normal",
                         type = "apeglm")

# convert to dataframe
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

# how many DEGs?
sig <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
message("upregulated: ", sum(sig$log2FoldChange > 0))
message("downregulated: ", sum(sig$log2FoldChange < 0))

# save results
write.csv(res_df, "results/DEG_tumor_vs_normal_full.csv", row.names = FALSE)
write.csv(sig,    "results/DEG_tumor_vs_normal_significant.csv", row.names = FALSE)


# ── Volcano plot ──────────────────────────────────────────────────────────────

res_plot <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    category = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# top genes to label - picking by smallest padj within sig
top_label <- res_plot %>%
  filter(category != "NS") %>%
  slice_min(padj, n = 20)

volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj),
                                colour = category)) +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_point(data = top_label, size = 2.2) +
  geom_text_repel(data = top_label, aes(label = gene),
                  size = 3, max.overlaps = 15,
                  box.padding = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey75")) +
  labs(
    title    = "Differential Expression: HCC Tumor vs Normal",
    subtitle = "FDR < 0.05  |  |log2FC| > 1",
    x        = "log2 Fold Change",
    y        = "-log10 adjusted p-value",
    colour   = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("plots/volcano_HCC.pdf", volcano, width = 8, height = 7)
ggsave("plots/volcano_HCC.png", volcano, width = 8, height = 7, dpi = 300)


# ── Heatmap of top DEGs ───────────────────────────────────────────────────────

top50 <- sig %>% slice_min(padj, n = 50) %>% pull(gene)

vst_data <- vst(dds, blind = FALSE)
mat      <- assay(vst_data)[top50, ]
mat      <- mat - rowMeans(mat)   # row-center

anno_col <- data.frame(
  condition = meta$condition,
  row.names = rownames(meta)
)

ann_colors <- list(condition = c(tumor = "#d62728", normal = "#1f77b4"))

pheatmap(mat,
         annotation_col  = anno_col,
         annotation_colors = ann_colors,
         show_rownames   = TRUE,
         fontsize_row    = 7,
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         color           = colorRampPalette(c("#053061","white","#67001f"))(100),
         border_color    = NA,
         main            = "Top 50 DEGs — HCC Tumor vs Normal (VST, row-centred)",
         filename        = "plots/heatmap_top50DEGs.pdf",
         width = 10, height = 12)


# ── PCA ───────────────────────────────────────────────────────────────────────

pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c(tumor = "#d62728", normal = "#1f77b4")) +
  labs(
    title  = "PCA — VST normalised counts",
    x      = paste0("PC1: ", pct_var[1], "% variance"),
    y      = paste0("PC2: ", pct_var[2], "% variance"),
    colour = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/PCA_HCC.pdf", pca_plot, width = 7, height = 6)
