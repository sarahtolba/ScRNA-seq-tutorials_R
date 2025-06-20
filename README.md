# 🔬 Comprehensive scRNA-seq Analysis in R using Seurat

This repository is an educational resource for researchers and students to learn **single-cell RNA-seq (scRNA-seq) analysis** using the **Seurat** package in R. It covers the entire standard pipeline and advanced topics, supporting multiple input formats including `.h5`, `.mtx`, and `.csv`.

---

## 🎯 Goals

- Teach single-cell RNA-seq analysis using R/Seurat
- Cover both **basic** and **advanced** workflows
- Serve as a reproducible and flexible template for various datasets

---

## 📚 Features Covered

✅ Basic Seurat workflow (normalization, clustering, visualization)  
✅ Support for input formats (`.h5`, `.mtx`)  
✅ Marker gene identification and cluster annotation  
✅ Cell type annotation using known markers or tools (e.g., SingleR)  
✅ Pseudotime analysis (e.g., with Monocle3 or Slingshot)  
✅ Cell–cell communication analysis (e.g., with CellChat)  
✅ Doublet detection and removal using DoubletFinder  
✅ Integration of multiple datasets  
✅ Clean and customizable plots for publications  

---

## 🛠️ Requirements

- R ≥ 4.2
- R packages:  
  `Seurat`, `patchwork`, `dplyr`, `ggplot2`,  
  `DoubletFinder`, `CellChat`, `monocle3` or `slingshot`, `SingleR`, etc.

To install core packages:

```r
install.packages("Seurat")
# Additional packages:
# BiocManager::install("monocle3") or install.packages("slingshot")
