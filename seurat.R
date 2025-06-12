# Load required libraries
library(Seurat)
library(tidyverse)

#---------------------------------------------------------------------
# Load data (10X Genomics .h5 format)
pbmcs_s_m = Read10X_h5(
  filename = "../Desktop/5k_Human_Donor3_PBMC_3p_gem-x_5k_Human_Donor3_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5"
)
str(pbmcs_s_m)
class(pbmcs_s_m)
dim(pbmcs_s_m)

#---------------------------------------------------------------------
# Let's make the Seurat object
seurat_object = CreateSeuratObject(
  counts = pbmcs_s_m,
  project = "PBMCS",
  min.cells = 3,
  min.features = 200
)
str(seurat_object)
view(seurat_object)

head(seurat_object@meta.data)

# % MT reads
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") 
# associated with low-quality reads/mitochondrial contamination
View(seurat_object@meta.data)
# percent.mt is the percentage of RNA transcripts (UMIs) in a cell that come from mitochondrial genes.

#---------------------------------------------------------------------
# Visualize quality control metrics
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
# to visualize the columns from meta data together

# A good-quality dataset should follow a straight line
FeatureScatter(seurat_object , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# because good quality cells not just have good numbers of genes but also counts

#---------------------------------------------------------------------
# 2. Filtering
seurat_object.filtered <- subset(
  seurat_object,
  subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5
)

VlnPlot(seurat_object.filtered , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_object.filtered , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') 

# upper left corner = genes present but not enough depth
# bottom right = droplets or outliers, maybe damaged cells

#---------------------------------------------------------------------
# 3. Normalize data
# default values: seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
seurat_object <- NormalizeData(seurat_object)
str(seurat_object) 
# all the commands after filtering will be in the command slot 

#---------------------------------------------------------------------
# 4. Identify highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#---------------------------------------------------------------------
# 5. Scaling
# to modulate technical noise: batch effect or biological differences (like cell cycle)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
# because I will perform linear dimensionality reduction 
# It ensures that highly expressed genes don’t dominate the analysis just because of their raw magnitude.
# It puts all genes on the same scale, so PCA, clustering, etc., are fair and meaningful.

str(seurat_object)

#---------------------------------------------------------------------
# 6. Perform linear dimensionality reduction
# to identify sources of heterogeneity in the data
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Visualize PCA results
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

# Determine dimensionality of the data
ElbowPlot(seurat_object) 
# These are the principal components and they are ranked according to the variance.
# I will pick the PCA where the dots start to settle down because there's little additional variance after that.

#---------------------------------------------------------------------
# 7. Clustering
# Gather the cells that have similar expression patterns
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)

# Understanding resolution
seurat_object <- FindClusters(seurat_object, resolution = c(0.1, 0.3, 0.5, 0.7, 1)) 
# The higher the resolution number, the more clusters you'll get.

View(seurat_object@meta.data)

# Check which resolution works the best to separate the clusters
DimPlot(seurat_object, group.by = "RNA_snn_res.0.1", label = TRUE)

# Set identity of clusters
Idents(seurat_object) <- "RNA_snn_res.0.1"

#---------------------------------------------------------------------
# 8. Non-linear dimensionality reduction
# Group cells of different dimensions together


seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Note: you can set `label = TRUE` or use the LabelClusters function to label individual clusters
DimPlot(seurat_object, reduction = "umap")

# don't change anything, just organize it
