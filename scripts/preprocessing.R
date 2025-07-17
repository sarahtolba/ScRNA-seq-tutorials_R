#---------------------------------------------------------------------
# 1. Load Required Libraries
#---------------------------------------------------------------------
library(Seurat)
library(tidyverse)

#---------------------------------------------------------------------
# 2. Load 10X Genomics Data (.h5 format)
#---------------------------------------------------------------------
pbmcs_s_m <- Read10X_h5(
  filename = "../data/project1/5k_Human_Donor3_PBMC_3p_gem-x_5k_Human_Donor3_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5"
)
dim(pbmcs_s_m)        # Check dimensions
class(pbmcs_s_m)      # Should be list or matrix

#---------------------------------------------------------------------
# 3. Create Seurat Object
#---------------------------------------------------------------------
seurat_object <- CreateSeuratObject(
  counts = pbmcs_s_m,
  project = "PBMCS",
  min.cells = 3,
  min.features = 200
)
head(seurat_object@meta.data)

#---------------------------------------------------------------------
# 4. Calculate % Mitochondrial Reads
#---------------------------------------------------------------------
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

#---------------------------------------------------------------------
# 5. QC Visualizations
#---------------------------------------------------------------------
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

#---------------------------------------------------------------------
# 6. Filter Low-quality Cells
#---------------------------------------------------------------------
seurat_object.filtered <- subset(
  seurat_object,
  subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5
)

# QC after filtering
VlnPlot(seurat_object.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_object.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

#---------------------------------------------------------------------
# 7. Normalize Data
#---------------------------------------------------------------------
seurat_object.filtered <- NormalizeData(seurat_object.filtered)

#---------------------------------------------------------------------
# 8. Identify Highly Variable Features
#---------------------------------------------------------------------
seurat_object.filtered <- FindVariableFeatures(
  seurat_object.filtered,
  selection.method = "vst",
  nfeatures = 2000
)

# Plot HVGs
top10 <- head(VariableFeatures(seurat_object.filtered), 10)
plot1 <- VariableFeaturePlot(seurat_object.filtered)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#---------------------------------------------------------------------
# 9. Scale Data
#---------------------------------------------------------------------
all.genes <- rownames(seurat_object.filtered)
seurat_object.filtered <- ScaleData(seurat_object.filtered, features = all.genes)

#---------------------------------------------------------------------
# 10. PCA (Linear Dimensionality Reduction)
#---------------------------------------------------------------------
seurat_object.filtered <- RunPCA(seurat_object.filtered, features = VariableFeatures(seurat_object.filtered))

print(seurat_object.filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_object.filtered, dims = 1, cells = 500, balanced = TRUE)

# Determine optimal number of PCs
ElbowPlot(seurat_object.filtered)

#---------------------------------------------------------------------
# 11. Clustering
#---------------------------------------------------------------------
seurat_object.filtered <- FindNeighbors(seurat_object.filtered, dims = 1:10)
seurat_object.filtered <- FindClusters(seurat_object.filtered, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# Visualize clustering at different resolutions
DimPlot(seurat_object.filtered, group.by = "RNA_snn_res.0.1", label = TRUE)

# Set default identity
Idents(seurat_object.filtered) <- "RNA_snn_res.0.1"

#---------------------------------------------------------------------
# 12. UMAP (Non-linear Dimensionality Reduction)
#---------------------------------------------------------------------
seurat_object.filtered <- RunUMAP(seurat_object.filtered, dims = 1:10)

# UMAP Plot
DimPlot(seurat_object.filtered, reduction = "umap", label = TRUE)

#save seurat object 
saveRDS(seurat_object.filtered, file = "../seurat objects/seurat_object_filtered.rds")

