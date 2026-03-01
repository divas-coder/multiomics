# pathway enrichment for HCC DEGs
# ORA using clusterProfiler - GO and KEGG
# also added GSEA at the bottom
# DM - updated march 2023 to add reactome

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(ReactomePA)

# load DEG results from previous step
# run DEG_analysis_HCC.R first if you havent already
deg <- read.csv("results/DEG_tumor_vs_normal_significant.csv")

# separate up and down regulated
up_genes   <- deg %>% filter(log2FoldChange > 0) %>% pull(gene)
down_genes <- deg %>% filter(log2FoldChange < 0) %>% pull(gene)

# background = all tested genes
all_deg <- read.csv("results/DEG_tumor_vs_normal_full.csv")
background <- all_deg$gene

message("upregulated genes: ", length(up_genes))
message("downregulated genes: ", length(down_genes))


# convert symbols to entrez IDs - needed for most databases
convert_to_entrez <- function(gene_symbols) {
  bitr(gene_symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db)
}

up_entrez  <- convert_to_entrez(up_genes)
dn_entrez  <- convert_to_entrez(down_genes)
bg_entrez  <- convert_to_entrez(background)

# how many converted?
message(round(nrow(up_entrez) / length(up_genes) * 100, 1), "% upregulated genes mapped")
message(round(nrow(dn_entrez) / length(down_genes) * 100, 1), "% downregulated genes mapped")


# ── GO enrichment ──────────────────────────────────────────────────────────────

run_go_ora <- function(entrez_ids, bg_entrez, ont = "BP", label = "") {
  enrichGO(
    gene          = entrez_ids$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}

go_up_bp <- run_go_ora(up_entrez, bg_entrez, "BP", "upregulated")
go_dn_bp <- run_go_ora(dn_entrez, bg_entrez, "BP", "downregulated")

# how many terms?
if (!is.null(go_up_bp)) message("GO BP upregulated: ", nrow(go_up_bp), " terms")
if (!is.null(go_dn_bp)) message("GO BP downregulated: ", nrow(go_dn_bp), " terms")

# dotplots
if (!is.null(go_up_bp) && nrow(go_up_bp) > 0) {
  p_up <- dotplot(go_up_bp, showCategory = 20) +
    labs(title = "GO Biological Process — Upregulated genes (HCC)") +
    theme(plot.title = element_text(face = "bold", size = 11))
  ggsave("plots/GO_BP_upregulated.pdf", p_up, width = 9, height = 8)
}

if (!is.null(go_dn_bp) && nrow(go_dn_bp) > 0) {
  p_dn <- dotplot(go_dn_bp, showCategory = 20) +
    labs(title = "GO Biological Process — Downregulated genes (HCC)") +
    theme(plot.title = element_text(face = "bold", size = 11))
  ggsave("plots/GO_BP_downregulated.pdf", p_dn, width = 9, height = 8)
}

# save tables
write.csv(as.data.frame(go_up_bp), "results/GO_BP_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(go_dn_bp), "results/GO_BP_downregulated.csv", row.names = FALSE)


# ── KEGG enrichment ─────────────────────────────────────────────────────────────

kegg_up <- enrichKEGG(
  gene          = up_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 10
)

kegg_dn <- enrichKEGG(
  gene          = dn_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 10
)

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  ggsave("plots/KEGG_upregulated.pdf",
         dotplot(kegg_up, showCategory = 20) +
           labs(title = "KEGG Pathways — Upregulated"),
         width = 9, height = 7)
  write.csv(as.data.frame(kegg_up), "results/KEGG_upregulated.csv", row.names = FALSE)
}

if (!is.null(kegg_dn) && nrow(kegg_dn) > 0) {
  ggsave("plots/KEGG_downregulated.pdf",
         dotplot(kegg_dn, showCategory = 20) +
           labs(title = "KEGG Pathways — Downregulated"),
         width = 9, height = 7)
  write.csv(as.data.frame(kegg_dn), "results/KEGG_downregulated.csv", row.names = FALSE)
}


# ── Reactome ───────────────────────────────────────────────────────────────────
# added this later - reactome often gives cleaner results than KEGG for cancer

reactome_up <- enrichPathway(
  gene          = up_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

if (!is.null(reactome_up) && nrow(reactome_up) > 0) {
  ggsave("plots/Reactome_upregulated.pdf",
         dotplot(reactome_up, showCategory = 20) +
           labs(title = "Reactome Pathways — Upregulated (HCC)"),
         width = 10, height = 8)
  write.csv(as.data.frame(reactome_up),
            "results/Reactome_upregulated.csv", row.names = FALSE)
}


# ── GSEA ───────────────────────────────────────────────────────────────────────
# using stat column for ranking - more stable than just lfc

gsea_input <- all_deg %>%
  filter(!is.na(stat)) %>%
  distinct(gene, .keep_all = TRUE)

# merge with entrez IDs
gsea_map <- convert_to_entrez(gsea_input$gene)
gsea_merged <- merge(gsea_input, gsea_map, by.x = "gene", by.y = "SYMBOL")

ranked_list <- setNames(gsea_merged$stat, gsea_merged$ENTREZID)
ranked_list <- sort(ranked_list, decreasing = TRUE)

# GSEA GO
gsea_go <- gseGO(
  geneList     = ranked_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  seed         = 2023
)

if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  message("GSEA GO enriched terms: ", nrow(gsea_go))

  # enrichment map
  gsea_go2 <- pairwise_termsim(gsea_go)
  p_emap <- emapplot(gsea_go2, showCategory = 30) +
    labs(title = "GSEA Enrichment Map — GO BP")
  ggsave("plots/GSEA_GO_enrichmap.pdf", p_emap, width = 12, height = 10)

  # top pathway running score plots
  top_ids <- head(gsea_go@result$ID, 5)
  for (pid in top_ids) {
    p_score <- gseaplot2(gsea_go, geneSetID = pid, title = pid)
    fname   <- paste0("plots/GSEA_running_score_", gsub(":", "_", pid), ".pdf")
    ggsave(fname, p_score, width = 9, height = 5)
  }

  write.csv(as.data.frame(gsea_go), "results/GSEA_GO_results.csv", row.names = FALSE)
}

# GSEA KEGG
gsea_kegg <- gseKEGG(
  geneList     = ranked_list,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  seed         = 2023
)

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  ggsave("plots/GSEA_KEGG_dotplot.pdf",
         dotplot(gsea_kegg, showCategory = 20, split = ".sign") +
           facet_grid(. ~ .sign) +
           labs(title = "GSEA KEGG — Activated vs Suppressed pathways"),
         width = 14, height = 8)
  write.csv(as.data.frame(gsea_kegg), "results/GSEA_KEGG_results.csv", row.names = FALSE)
}

message("pathway analysis done - check plots/ and results/ folders")
