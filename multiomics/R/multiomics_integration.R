# multiomics integration - combining transcriptomics + microbiome
# using MOFA2 for factor analysis
# also doing correlation analysis between microbial abundance and gene expression
# DM - this is the final integration step, run after individual analyses

library(MOFA2)
library(ggplot2)
library(dplyr)
library(corrplot)
library(Hmisc)   # for rcorr function

# load results from individual analyses
deg_res  <- read.csv("results/DEG_tumor_vs_normal_full.csv")
mb_res   <- read.csv("results/microbiome_differential_abundance.csv")

# load normalised expression matrix and microbiome abundance
expr_norm <- read.csv("data/TCGA_LIHC_normalised_expr.csv", row.names = 1)
mb_abund  <- read.csv("data/microbiome_genus_abundance.csv", row.names = 1)

# align samples
shared_samples <- intersect(colnames(expr_norm), colnames(mb_abund))
expr_shared    <- expr_norm[, shared_samples]
mb_shared      <- mb_abund[, shared_samples]

message("shared samples for integration: ", length(shared_samples))


# ── MOFA2 factor analysis ──────────────────────────────────────────────────────
# using top variable features from each omics layer

# top 2000 most variable genes
var_genes <- apply(expr_shared, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:2000]

# top 100 most variable microbes
var_mb <- apply(log1p(mb_shared), 1, var)
top_mb <- names(sort(var_mb, decreasing = TRUE))[1:100]

# prepare MOFA input - list of matrices
# each matrix: features x samples
mofa_input <- list(
  transcriptomics = as.matrix(expr_shared[top_genes, ]),
  microbiome      = as.matrix(log1p(mb_shared[top_mb, ]))
)

# scale each view
mofa_input <- lapply(mofa_input, function(mat) {
  t(scale(t(mat)))  # scale features
})

# create MOFA object
mofa <- create_mofa(mofa_input)

# check data overview
plot_data_overview(mofa)

# set MOFA options
mofa_opts <- get_default_model_options(mofa)
mofa_opts$num_factors <- 10   # start with 10, can increase

train_opts <- get_default_training_options(mofa)
train_opts$seed           <- 42
train_opts$convergence_mode <- "fast"

mofa <- prepare_mofa(mofa,
                      model_options    = mofa_opts,
                      training_options = train_opts)

# train model - takes a few minutes
mofa_trained <- run_mofa(mofa, use_basilisk = TRUE)

# save trained model
saveRDS(mofa_trained, "results/MOFA2_trained_model.rds")

# load sample metadata
meta <- read.csv("data/TCGA_LIHC_clinical.csv", row.names = 1)
meta_shared <- meta[shared_samples, , drop = FALSE]
samples_metadata(mofa_trained) <- meta_shared

# plot variance explained per factor
p_var <- plot_variance_explained(mofa_trained, max.var = 15) +
  labs(title = "MOFA2 — Variance Explained per Factor")
ggsave("plots/MOFA2_variance_explained.pdf", p_var, width = 8, height = 5)

# plot factors coloured by condition
p_factors <- plot_factor(mofa_trained, factors = 1:4,
                          color_by = "condition") +
  scale_colour_manual(values = c(tumor = "#d62728", normal = "#1f77b4"))
ggsave("plots/MOFA2_factors.pdf", p_factors, width = 10, height = 8)

# top feature weights for Factor 1
p_weights_rna <- plot_top_weights(mofa_trained,
                                   view   = "transcriptomics",
                                   factor = 1,
                                   nfeatures = 25) +
  labs(title = "Top transcriptomic features — Factor 1")

p_weights_mb <- plot_top_weights(mofa_trained,
                                  view   = "microbiome",
                                  factor = 1,
                                  nfeatures = 25) +
  labs(title = "Top microbiome features — Factor 1")

ggsave("plots/MOFA2_weights_transcriptomics_F1.pdf", p_weights_rna, width = 7, height = 8)
ggsave("plots/MOFA2_weights_microbiome_F1.pdf",      p_weights_mb,  width = 7, height = 8)


# ── Transcriptome-Microbiome Correlation ────────────────────────────────────────
# checking if specific genera correlate with DEG expression
# focusing on sig DEGs and sig DA genera

sig_deg_genes <- deg_res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(gene)

sig_da_genera <- mb_res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(Genus) %>%
  na.omit()

# subset matrices to sig features
expr_sub <- expr_shared[sig_deg_genes[sig_deg_genes %in% rownames(expr_shared)], ]
mb_sub   <- log1p(mb_shared[sig_da_genera[sig_da_genera %in% rownames(mb_shared)], ])

if (nrow(expr_sub) > 0 && nrow(mb_sub) > 0) {

  # spearman correlations - more robust than pearson for this
  cor_mat <- cor(t(expr_sub), t(mb_sub), method = "spearman")

  # p-values using rcorr
  cor_pmat <- rcorr(t(as.matrix(expr_sub)), t(as.matrix(mb_sub)),
                     type = "spearman")

  # keep only gene-microbe pairs (not gene-gene or microbe-microbe)
  n_genes  <- nrow(expr_sub)
  n_genera <- nrow(mb_sub)
  cor_pairs <- cor_pmat$r[1:n_genes, (n_genes+1):(n_genes+n_genera)]
  pval_pairs <- cor_pmat$P[1:n_genes, (n_genes+1):(n_genes+n_genera)]

  # FDR correction
  pval_adj <- matrix(p.adjust(as.vector(pval_pairs), method = "BH"),
                      nrow = nrow(pval_pairs))
  rownames(pval_adj) <- rownames(cor_pairs)
  colnames(pval_adj) <- colnames(cor_pairs)

  # significant pairs
  sig_mask <- pval_adj < 0.05 & abs(cor_pairs) > 0.3

  message("significant gene-microbe correlations: ", sum(sig_mask, na.rm = TRUE))

  # heatmap of top correlated pairs
  sig_genes_any  <- rownames(cor_pairs)[rowSums(sig_mask, na.rm = TRUE) > 0]
  sig_genera_any <- colnames(cor_pairs)[colSums(sig_mask, na.rm = TRUE) > 0]

  if (length(sig_genes_any) > 0 && length(sig_genera_any) > 0) {
    # take top 30 genes by number of significant correlations
    top_genes_corr  <- sig_genes_any[order(rowSums(sig_mask[sig_genes_any, , drop = FALSE],
                                                     na.rm = TRUE),
                                            decreasing = TRUE)][1:min(30, length(sig_genes_any))]

    cor_sub_plot <- cor_pairs[top_genes_corr, sig_genera_any, drop = FALSE]

    pdf("plots/transcriptome_microbiome_correlation_heatmap.pdf", width = 12, height = 10)
    corrplot(cor_sub_plot,
             method     = "color",
             col        = colorRampPalette(c("#053061", "white", "#67001f"))(200),
             tl.cex     = 0.7,
             tl.col     = "black",
             cl.lim     = c(-1, 1),
             addCoef.col = NULL,
             title      = "Spearman Correlations: DEGs vs DA Genera",
             mar        = c(0, 0, 2, 0))
    dev.off()

    # save correlation results
    cor_long <- as.data.frame(as.table(cor_pairs)) %>%
      rename(gene = Var1, genus = Var2, spearman_r = Freq) %>%
      mutate(
        padj = as.vector(pval_adj),
        significant = padj < 0.05 & abs(spearman_r) > 0.3
      ) %>%
      filter(significant) %>%
      arrange(padj)

    write.csv(cor_long, "results/transcriptome_microbiome_correlations.csv",
              row.names = FALSE)

    message("saved ", nrow(cor_long), " significant gene-microbe pairs")
  }
}

message("multiomics integration complete")
message("check plots/ for MOFA2 and correlation visualisations")
