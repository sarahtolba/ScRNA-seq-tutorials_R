
# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)
#determine data location
loc = list.dirs(path = "../Desktop/Rprojects/scRNA-seq/project/data/", full.names = FALSE,recursive = FALSE)
base_path = "../Desktop/Rprojects/scRNA-seq/project/data/"
#we have 4 folders for for sampples so we will
#creat a loop to do all steps for each folder at one time


 for (i in loc){
   
   sample_name <- i
   
   cts = ReadMtx(mtx = paste0(base_path,i,"/matrix.mtx"),
           features = paste0(base_path,i,"/features.tsv"),
           cells =paste0(base_path,i,"/barcodes.tsv"))
   
   # create seurat objects
   assign(sample_name, CreateSeuratObject(counts = cts))
   
  
 }  


# merin and not interation so we can berform the standard work flow
merged_seurat <- merge(
  x = control,
  y = c(phenoformin, polyIC, combination),
  add.cell.ids = c("control", "phenoformin", "polyIC", "combination"),
  project = "CII"
)

#qc
view(merged_seurat@meta.data)

merged_seurat@meta.data$sample = row.names(merged_seurat@meta.data)
merged_seurat@meta.data = separate(merged_seurat@meta.data, col = 'sample', into = c('Type', 'Barcode'), 
                                   sep = '_')
#to ggmake sure that our data present
unique(merged_seurat@meta.data$Type)


# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
unique(merged_seurat@meta.data$mitoPercent)

# Visualize quality control metrics
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3) 

# high mitochondrial content cells had already been removed prior to analysis

# A good-quality dataset should follow a straight line
 FeatureScatter(merged_seurat , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

#filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)
# needs droplet finders 
merged_seurat_filtered

merged_seurat



 FeatureScatter(merged_seurat_filtered , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

#in order to perform interation we need to see for batch effect in a low dimention 
#so we need to perform the standard work flow


merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

# Visualize PCA results
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

 p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type')




 
 
 
 
 
 
