# data : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272922
# ===============================
# 📦 Load Required Libraries
# ===============================
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)

# ===============================
# 📁 Set Paths and Sample Info
# ===============================
base_path <- "../Desktop/Rprojects/scRNA-seq/project/data/"
sample_dirs <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)

# ===============================
# 🔁 Read Each Sample and Create Seurat Object
# ===============================
for (sample in sample_dirs) {
  counts <- ReadMtx(
    mtx = paste0(base_path, sample, "/matrix.mtx"),
    features = paste0(base_path, sample, "/features.tsv"),
    cells = paste0(base_path, sample, "/barcodes.tsv")
  )
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  assign(sample, seurat_obj)
}

# ===============================
# 🔗 Merge Samples (Without Integration)
# ===============================
merged_seurat <- merge(
  x = control,
  y = c(phenoformin, polyIC, combination),
  add.cell.ids = c("control", "phenoformin", "polyIC", "combination"),
  project = "CII"
)

# ===============================
# 🧬 Add Sample Type Metadata
# ===============================
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(
  merged_seurat@meta.data,
  col = "sample",
  into = c("Type", "Barcode"),
  sep = "_"
)

# Confirm sample types
unique(merged_seurat$Type)

# ===============================
# 🔍 Quality Control
# ===============================
# Calculate % mitochondrial gene expression
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)

# Feature scatter for count vs. features
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# ===============================
# 🧼 Filtering Low-Quality Cells
# ===============================
merged_seurat_filtered <- subset(
  merged_seurat,
  subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10
)

# Optional: Save filtered object for inspection
merged_seurat_filtered

# Feature scatter post-filtering
FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# ===============================
# 🧪 Standard Seurat Workflow (Before Integration)
# ===============================
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)

# PCA inspection
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(merged_seurat_filtered)

# ===============================
# 🔗 Clustering & UMAP
# ===============================
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)

# UMAP colored by sample type
DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Type")
