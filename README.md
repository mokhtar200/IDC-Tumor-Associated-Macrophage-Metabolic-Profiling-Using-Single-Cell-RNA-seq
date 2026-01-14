# Tumor-Associated Macrophage (TAM) Profiling in IDC Breast Cancer scRNA-seq

This repository contains a pipeline for analyzing **Tumor-Associated Macrophages (TAMs)** from **IDC breast cancer single-cell RNA-seq data** using Seurat and SingleR.

---

## Features

- Quality control and SCTransform normalization
- Dimensionality reduction (PCA, UMAP)
- Automatic cell-type annotation using **SingleR** with **HumanPrimaryCellAtlasData**
- Extraction of **immune cells** and **TAMs** (Macrophages + Monocytes)
- Functional marker visualization (M1, M2, APC markers)
- Metabolic scoring using **FDCA (Flux Descriptive Cellular Automata)**
- Safe handling of small TAM subsets and missing genes
- Results saved as CSV summary tables

---

## Requirements

- R >= 4.2
- Packages: `Seurat`, `tidyverse`, `patchwork`, `SingleR`, `celldex`, `scales`

Install packages if needed:

```r
install.packages(c("Seurat","tidyverse","patchwork","scales"))
BiocManager::install(c("SingleR","celldex"))
