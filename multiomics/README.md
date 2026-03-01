# 🧬 Multiomics Analysis — HCC

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-Differential_Expression-E74C3C?style=flat-square)](https://bioconductor.org/packages/DESeq2)
[![MOFA2](https://img.shields.io/badge/MOFA2-Multi--Omics_Integration-8E44AD?style=flat-square)](https://biofam.github.io/MOFA2/)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)

Integrative transcriptomic and microbiome analysis in **Hepatocellular Carcinoma (HCC)** — identifying prognostic biomarkers, dysregulated pathways, and gut microbiome alterations associated with liver cancer.

This work supports the findings in:

> Mishra D., Mishra A., Rai S.N., Vamanu E., Singh M.P. *"Identification of prognostic biomarkers for suppressing tumorigenesis and metastasis of hepatocellular carcinoma through transcriptome analysis."* **Diagnostics**, 2023. https://doi.org/10.3390/diagnostics13050965

---

## 📋 Analysis Modules

| Script | Description |
|--------|-------------|
| `R/DEG_analysis_HCC.R` | Differential gene expression using DESeq2 |
| `R/pathway_enrichment_HCC.R` | ORA (GO/KEGG/Reactome) and GSEA |
| `R/survival_analysis_HCC.R` | Kaplan-Meier and Cox PH modelling (TCGA-LIHC) |
| `R/microbiome_analysis_HCC.R` | 16S microbiome diversity and differential abundance |
| `R/multiomics_integration.R` | MOFA2 factor analysis + transcriptome-microbiome correlations |

---

## 🗂️ Project Structure

```
multiomics/
├── R/
│   ├── DEG_analysis_HCC.R
│   ├── pathway_enrichment_HCC.R
│   ├── survival_analysis_HCC.R
│   ├── microbiome_analysis_HCC.R
│   └── multiomics_integration.R
├── data/                     # input data (not tracked in git)
├── plots/                    # output figures
├── results/                  # output tables
└── docs/
```

---

## ⚙️ Dependencies

```r
# Bioconductor
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot",
                        "ReactomePA", "org.Hs.eg.db", "phyloseq",
                        "microbiome", "MOFA2"))

# CRAN
install.packages(c("survival", "survminer", "vegan",
                    "corrplot", "Hmisc", "ggrepel",
                    "ggplot2", "dplyr", "tidyr"))
```

---

## 🚀 How to Run

Run scripts in order:

```r
source("R/DEG_analysis_HCC.R")           # Step 1 — DEG
source("R/pathway_enrichment_HCC.R")     # Step 2 — Pathways
source("R/survival_analysis_HCC.R")      # Step 3 — Survival
source("R/microbiome_analysis_HCC.R")    # Step 4 — Microbiome
source("R/multiomics_integration.R")     # Step 5 — Integration
```

---

## 📊 Key Outputs

- Volcano plots, MA plots, DEG heatmaps
- GO/KEGG/Reactome dotplots and GSEA enrichment maps
- Kaplan-Meier survival curves + Cox forest plots
- Bray-Curtis PCoA, alpha diversity boxplots
- MOFA2 factor plots + transcriptome-microbiome correlation heatmap

---

## 👤 Author

**Dr. Divya Mishra** — Ph.D. Bioinformatics

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=flat-square&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/dr-divya-mishra/)
[![ResearchGate](https://img.shields.io/badge/ResearchGate-00CCBB?style=flat-square&logo=researchgate&logoColor=white)](https://www.researchgate.net/profile/Divya-Mishra-13)
[![Email](https://img.shields.io/badge/Email-D14836?style=flat-square&logo=gmail&logoColor=white)](mailto:mishra.divya76@gmail.com)
