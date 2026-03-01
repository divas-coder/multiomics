# microbiome analysis - 16S rRNA gut microbiome data
# part of the multiomics integration work
# comparing gut microbiome between HCC patients and healthy controls
# DM - this was tricky to get right, took a few iterations

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(microbiome)  # handy wrappers around phyloseq
library(DESeq2)      # using for differential abundance too

# load phyloseq object
# created from QIIME2 output - see docs/QIIME2_processing_notes.txt
ps <- readRDS("data/phyloseq_HCC_microbiome.rds")

# quick look
ps
sample_data(ps)$condition %>% table()

# ── Alpha diversity ────────────────────────────────────────────────────────────

alpha_div <- microbiome::alpha(ps, index = c("shannon", "observed", "chao1"))
alpha_div$condition <- sample_data(ps)$condition

# test differences - using wilcoxon since alpha div rarely normally distributed
shannon_test <- wilcox.test(shannon ~ condition, data = alpha_div)
observed_test <- wilcox.test(observed ~ condition, data = alpha_div)

message("Shannon p-value: ", round(shannon_test$p.value, 4))
message("Observed richness p-value: ", round(observed_test$p.value, 4))

# plot alpha diversity
alpha_long <- alpha_div %>%
  tidyr::pivot_longer(cols = c(shannon, observed, chao1),
                       names_to = "metric", values_to = "value")

p_alpha <- ggplot(alpha_long, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c(HCC = "#d62728", control = "#1f77b4")) +
  labs(
    title = "Alpha diversity — HCC vs Control",
    x = NULL, y = "Diversity index"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )

ggsave("plots/alpha_diversity.pdf", p_alpha, width = 9, height = 5)


# ── Beta diversity ─────────────────────────────────────────────────────────────

# rarefy to even depth first - use min library size
min_depth <- min(sample_sums(ps))
message("rarefying to: ", min_depth, " reads")
ps_rare <- rarefy_even_depth(ps, sample.size = min_depth, rngseed = 42)

# Bray-Curtis dissimilarity
bc_dist <- phyloseq::distance(ps_rare, method = "bray")

# PERMANOVA
pm_test <- adonis2(bc_dist ~ condition,
                    data   = as(sample_data(ps_rare), "data.frame"),
                    permutations = 999)
print(pm_test)

# Betadisper - test for homogeneity of dispersion
bd <- betadisper(bc_dist,
                  group = sample_data(ps_rare)$condition)
bd_test <- permutest(bd, permutations = 999)
print(bd_test)

# PCoA plot
pcoa <- ordinate(ps_rare, method = "PCoA", distance = "bray")
pcoa_df <- as.data.frame(pcoa$vectors[, 1:2])
colnames(pcoa_df) <- c("Axis1", "Axis2")
pcoa_df$condition <- sample_data(ps_rare)$condition
pcoa_df$sample    <- rownames(pcoa_df)

eig <- pcoa$values$Eigenvalues
pct <- round(eig / sum(eig) * 100, 1)

p_pcoa <- ggplot(pcoa_df, aes(Axis1, Axis2, colour = condition)) +
  geom_point(size = 3.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = c(HCC = "#d62728", control = "#1f77b4")) +
  labs(
    title    = "Beta diversity — Bray-Curtis PCoA",
    subtitle = paste0("PERMANOVA R²=", round(pm_test$R2[1], 3),
                      "  p=", pm_test$`Pr(>F)`[1]),
    x        = paste0("PCoA1 (", pct[1], "%)"),
    y        = paste0("PCoA2 (", pct[2], "%)"),
    colour   = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/beta_diversity_PCoA.pdf", p_pcoa, width = 7, height = 6)


# ── Differential abundance ─────────────────────────────────────────────────────
# using DESeq2 on genus level - more power than ANCOM for this dataset size

# agglomerate to genus level
ps_genus <- tax_glom(ps, taxrank = "Genus")

# remove genera with too many zeros - at least 20% prevalence
ps_genus_filt <- filter_taxa(ps_genus,
                               function(x) sum(x > 0) >= 0.2 * length(x),
                               TRUE)
message("genera after filtering: ", ntaxa(ps_genus_filt))

# convert to DESeq2
dds_mb <- phyloseq_to_deseq2(ps_genus_filt, ~ condition)
dds_mb <- DESeq(dds_mb, fitType = "local")

res_mb <- results(dds_mb,
                   contrast      = c("condition", "HCC", "control"),
                   alpha         = 0.05,
                   pAdjustMethod = "BH")

res_mb_df <- as.data.frame(res_mb) %>%
  tibble::rownames_to_column("ASV") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# add taxonomy
tax_df <- as.data.frame(tax_table(ps_genus_filt)) %>%
  tibble::rownames_to_column("ASV")

res_mb_annotated <- merge(res_mb_df, tax_df, by = "ASV") %>%
  arrange(padj)

sig_mb <- res_mb_annotated %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
message("differentially abundant genera: ", nrow(sig_mb))

write.csv(res_mb_annotated, "results/microbiome_differential_abundance.csv",
          row.names = FALSE)

# plot top DA genera
if (nrow(sig_mb) > 0) {
  top_genera <- sig_mb %>%
    slice_min(padj, n = min(30, nrow(sig_mb))) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "Enriched in HCC", "Depleted in HCC"),
           Genus = ifelse(is.na(Genus), paste0("Unknown (", Family, ")"), Genus))

  p_da <- ggplot(top_genera, aes(x = log2FoldChange,
                                   y = reorder(Genus, log2FoldChange),
                                   fill = direction)) +
    geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Enriched in HCC" = "#d62728",
                                  "Depleted in HCC" = "#1f77b4")) +
    labs(
      title = "Differentially Abundant Genera — HCC vs Control",
      x     = "log2 Fold Change",
      y     = "Genus",
      fill  = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave("plots/microbiome_DA_genera.pdf", p_da, width = 9, height = 8)
}

message("microbiome analysis done")
