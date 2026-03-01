# survival analysis for HCC - using TCGA-LIHC clinical data
# looking at whether top DEG expression correlates with patient survival
# DM - added this after reviewing the paper - reviewers asked for it

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(tidyr)

# load clinical data
# downloaded from TCGA data portal - TCGA-LIHC project
clinical <- read.csv("data/TCGA_LIHC_clinical.csv", row.names = 1)
expr     <- read.csv("data/TCGA_LIHC_normalised_expr.csv", row.names = 1)

# check dimensions
dim(clinical)
dim(expr)

# align samples
shared <- intersect(rownames(clinical), colnames(expr))
clinical <- clinical[shared, ]
expr     <- expr[, shared]

message(length(shared), " samples with both clinical and expression data")

# survival variables - OS and PFI
# OS = overall survival
# PFI = progression free interval
# both available in TCGA

# quick check on survival data
table(clinical$vital_status)
summary(clinical$OS.time)  # in days


# ── Kaplan-Meier for individual genes ──────────────────────────────────────────
# split into high/low by median expression

km_for_gene <- function(gene_name, clinical_df, expr_mat) {

  if (!gene_name %in% rownames(expr_mat)) {
    warning(gene_name, " not found in expression matrix")
    return(NULL)
  }

  gene_expr <- expr_mat[gene_name, ]

  df <- data.frame(
    sample    = names(gene_expr),
    expr      = as.numeric(gene_expr),
    OS_time   = clinical_df$OS.time / 30,   # convert days to months
    OS_status = ifelse(clinical_df$vital_status == "Dead", 1, 0)
  ) %>%
    filter(!is.na(OS_time), !is.na(OS_status), OS_time > 0)

  # median split
  med <- median(df$expr, na.rm = TRUE)
  df$group <- ifelse(df$expr >= med, "High", "Low")

  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = df)

  # log rank test
  lr  <- survdiff(Surv(OS_time, OS_status) ~ group, data = df)
  pval <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

  p <- ggsurvplot(
    fit,
    data         = df,
    pval         = TRUE,
    pval.method  = TRUE,
    conf.int     = TRUE,
    risk.table   = TRUE,
    palette      = c("#d62728", "#1f77b4"),
    legend.labs  = c(paste0("High (n=", sum(df$group == "High"), ")"),
                     paste0("Low (n=",  sum(df$group == "Low"),  ")")),
    title        = paste0(gene_name, " — OS in TCGA-LIHC"),
    xlab         = "Time (months)",
    ylab         = "Overall Survival probability",
    ggtheme      = theme_classic(base_size = 12)
  )

  list(plot = p, pval = pval, gene = gene_name)
}


# run KM for top prognostic genes from DEG analysis
# picked genes with literature support + strong FC
prognostic_genes <- c("TP53", "MKI67", "CCNB1", "CDK1", "AFP",
                       "GPC3", "EPCAM", "CD44", "ALDH1A1", "SOX9")

km_results <- list()
for (gene in prognostic_genes) {
  result <- km_for_gene(gene, clinical, expr)
  if (!is.null(result)) {
    km_results[[gene]] <- result
    fname <- paste0("plots/KM_", gene, "_OS.pdf")
    pdf(fname, width = 8, height = 7)
    print(result$plot)
    dev.off()
  }
}

# summary table of p-values
km_summary <- data.frame(
  gene = names(km_results),
  pval = sapply(km_results, function(x) x$pval)
) %>%
  arrange(pval) %>%
  mutate(significant = pval < 0.05)

print(km_summary)
write.csv(km_summary, "results/survival_pvalues_summary.csv", row.names = FALSE)


# ── Cox proportional hazards ────────────────────────────────────────────────────
# multivariate model adjusting for age, stage, grade

sig_genes_surv <- km_summary %>%
  filter(significant) %>%
  pull(gene)

if (length(sig_genes_surv) > 0) {
  message("Building Cox model with ", length(sig_genes_surv), " significant genes")

  # build expression matrix for Cox input
  expr_cox <- t(expr[sig_genes_surv, , drop = FALSE]) %>%
    as.data.frame()

  cox_df <- cbind(
    clinical[, c("OS.time", "vital_status", "age_at_index",
                 "ajcc_pathologic_stage", "tumor_grade")],
    expr_cox
  ) %>%
    mutate(
      OS_time   = OS.time / 30,
      OS_status = ifelse(vital_status == "Dead", 1, 0)
    ) %>%
    filter(!is.na(OS_time), !is.na(OS_status))

  # formula
  covars  <- c("age_at_index", "ajcc_pathologic_stage", sig_genes_surv)
  formula <- as.formula(paste("Surv(OS_time, OS_status) ~",
                               paste(covars, collapse = " + ")))

  cox_fit <- coxph(formula, data = cox_df)
  print(summary(cox_fit))

  # forest plot
  p_forest <- ggforest(cox_fit,
                        data    = cox_df,
                        main    = "Cox PH Model — OS in TCGA-LIHC",
                        fontsize = 0.85)
  ggsave("plots/cox_forest_plot.pdf", p_forest, width = 10, height = 8)

  # save cox results
  cox_res <- as.data.frame(summary(cox_fit)$coefficients)
  cox_res$gene <- rownames(cox_res)
  write.csv(cox_res, "results/cox_ph_results.csv", row.names = FALSE)
}

message("survival analysis complete")

# sessionInfo() saved just in case
# R version 4.3.1, survival 3.5-7, survminer 0.4.9
