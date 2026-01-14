suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SingleR)
  library(celldex)
  library(scales)
})

set.seed(123)

#Load Cell Ranger Data
#======================

data_dir <- "D:/Breast_Cancer_Single_Cell/filtered_feature_bc_matrix"

raw_counts <- Read10X(data.dir = data_dir)

sce <- CreateSeuratObject(
  counts = raw_counts,
  project = "IDC_scRNA",
  min.cells = 3,
  min.features = 200
)

#========================
# Quality Control
#========================

sce[["percent.mt"]] <- PercentageFeatureSet(
  sce,
  pattern = "^MT-"
)

sce <- subset(
  sce,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

#=============================
# Normalization (SCTransform)
#=============================

sce <- SCTransform(
  sce,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

DefaultAssay(sce) <- "SCT"

#==============================
# Dimensionality Reduction
#============================

sce <- RunPCA(sce, verbose = FALSE)

sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.4)

sce <- RunUMAP(sce, dims = 1:20)

DimPlot(sce, reduction = "umap", label = TRUE) + NoLegend()


#============================================
# Automatic Cell-Type Annotation (SingleR)
#==========================================

# Convert Seurat to SingleCellExperiment using SCT
sce_sce <- as.SingleCellExperiment(sce, assay = "SCT")

# Load reference
ref <- HumanPrimaryCellAtlasData()

# Run SingleR
pred <- SingleR(
  test = sce_sce,
  ref = ref,
  labels = ref$label.main
)

# Add labels back to Seurat
sce$SingleR_label <- pred$labels


#==========================
# Immune Cell Extraction (Observed Labels)
#=========================

immune <- subset(
  sce,
  subset = SingleR_label %in% c(
    "Macrophage",
    "Monocyte",
    "NK_cell",
    "DC"
  )
)

table(immune$SingleR_label)

DimPlot(
  immune,
  reduction = "umap",
  group.by = "SingleR_label",
  label = TRUE
)

#=================================
# Tumor-Associated Macrophage (TAM) Definition
#========================================
tam <- subset(
  immune,
  subset = SingleR_label %in% c("Macrophage", "Monocyte")
)

tam$cell_type <- "Tumor_Associated_Macrophage"

ncol(tam)


#=======================================
# Visual Sanity Check
#======================================
DimPlot(
  tam,
  reduction = "umap",
  group.by = "cell_type",
  pt.size = 2
) + NoLegend()


#=======================================
# TAM Functional Marker Panels
#=======================================
m1_markers <- c("IL1B", "TNF", "NOS2", "CXCL10")

m2_markers <- c("MRC1", "CD163", "ARG1", "IL10", "TGFB1")

apc_markers <- c("HLA-DRA", "HLA-DRB1", "CD74")

#==========================================
# TAM-Relevant Metabolic Programs
#==========================================
tam_metabolic_programs <- list(
  Glycolysis = c("HK2","PFKP","ALDOA","ENO1","PKM","LDHA"),
  OxPhos = c("NDUFA1","NDUFB8","SDHB","COX5A","ATP5F1A"),
  Fatty_Acid_Oxidation = c("CPT1A","ACADM","ACADVL"),
  Arginine_Metabolism = c("ARG1","NOS2","ASS1"),
  Tryptophan_Metabolism = c("IDO1"),
  Glutamine_Metabolism = c("GLS","SLC1A5")
)

tam_metabolic_programs <- lapply(
  tam_metabolic_programs,
  function(g) intersect(g, rownames(tam))
)

tam_metabolic_programs <- tam_metabolic_programs[
  sapply(tam_metabolic_programs, length) > 2
]
#==========================================
# FDCA: Descriptive Metabolic Scoring
#==========================================
tam <- AddModuleScore(
  tam,
  features = tam_metabolic_programs,
  assay = "SCT",
  name = "FDCA",
  nbin = 1   # For very small datasets, use 1 or 2
)


fdca_cols <- grep("^FDCA", colnames(tam@meta.data), value = TRUE)
names(fdca_cols) <- names(tam_metabolic_programs)

for (i in seq_along(fdca_cols)) {
  colnames(tam@meta.data)[
    colnames(tam@meta.data) == fdca_cols[i]
  ] <- names(tam_metabolic_programs)[i]
}


#==============================================
# Core Results: TAM Metabolic Landscape
#====================================================
FeaturePlot(
  tam,
  features = names(tam_metabolic_programs),
  reduction = "umap",
  ncol = 3,
  pt.size = 2
)

#==================================================
# Functional Marker Visualization
#==================================================
FeaturePlot(
  tam,
  features = c(m1_markers, m2_markers, apc_markers),
  reduction = "umap",
  ncol = 3,
  pt.size = 2
)


#===================================================
# SAFE Summary Table (No Statistics)
#================================================
tam_summary <- tam@meta.data %>%
  summarise(across(all_of(names(tam_metabolic_programs)), mean))

dir.create("results", showWarnings = FALSE)

write.csv(
  tam_summary,
  file = "results/TAM_Metabolic_Profile.csv",
  row.names = FALSE
)

