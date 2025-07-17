# Load required libraries
library(Seurat)
library(DoubletFinder)
library(ggplot2)

Reductions(seurat_obj)

# ------------------------------------------------------------------------------
# Step 1: Load pre-processed Seurat object
seurat_obj <- readRDS("../seurat objects/seurat_object_filtered.rds")

# ------------------------------------------------------------------------------
# Step 2: Run parameter sweep to find optimal pK
sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select best pK
best.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

# Estimate expected number of doublets
nExp_poi <- round(0.075 * ncol(seurat_obj))

# ------------------------------------------------------------------------------
# Step 3: Run DoubletFinder
seurat_obj <- doubletFinder(seurat_obj,
                            PCs = 1:20,
                            pN = 0.25,
                            pK = best.pK,
                            nExp = nExp_poi,
                            sct = FALSE)

# ------------------------------------------------------------------------------
# Step 4: Assign and visualize doublet status
df_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
seurat_obj$doublet_status <- seurat_obj@meta.data[[df_col]]

# UMAP plot colored by doublet status
DimPlot(seurat_obj, group.by = "doublet_status", reduction = "umap") +
  ggtitle("DoubletFinder: Singlets vs Doublets") +
  theme_minimal()

# ------------------------------------------------------------------------------
# Step 5: Filter to retain only singlets
seurat_obj_filtered <- subset(seurat_obj, subset = doublet_status == "Singlet")

# UMAP plot colored by doublet status
DimPlot(seurat_obj_filtered, group.by = "doublet_status", reduction = "umap") +
  ggtitle("DoubletFinder: Singlets vs Doublets") +
  theme_minimal()
# ------------------------------------------------------------------------------
# Step 6: Save the filtered object
saveRDS(seurat_obj_filtered, 
        file = "../seurat objects/seurat_obj_filtered.rds")
