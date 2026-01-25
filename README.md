# scClustEval <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Zaoqu-Liu/scClustEval/workflows/R-CMD-check/badge.svg)](https://github.com/Zaoqu-Liu/scClustEval/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/scClustEval)](https://CRAN.R-project.org/package=scClustEval)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/scClustEval)](https://cran.r-project.org/package=scClustEval)
<!-- badges: end -->

## Single Cell Clustering Evaluation and Optimization Framework

**scClustEval** is a comprehensive R package for evaluating and optimizing single-cell RNA-seq clustering results using self-projection machine learning approaches.

This is an R implementation inspired by the [SCCAF](https://github.com/SCCAF/sccaf) Python package.

## 🚀 Features

- **Clustering Assessment**: Evaluate clustering quality using self-projection with multiple ML classifiers
- **Clustering Optimization**: Iteratively merge poorly discriminated clusters
- **Multiple Classifiers**: Support for Logistic Regression, Random Forest, SVM, Naive Bayes, XGBoost, and more
- **Visualization**: ROC curves, confusion matrix heatmaps, Sankey diagrams
- **Seurat Integration**: Seamless workflow with Seurat objects (V4 & V5 compatible)
- **High Performance**: C++ acceleration for critical computations
- **Cross-Platform**: Works on Windows, macOS, and Linux

## 📦 Installation

### From GitHub (recommended)

```r
# Install devtools if not installed
if (!require("devtools")) install.packages("devtools")

# Install scClustEval
devtools::install_github("Zaoqu-Liu/scClustEval")
```

### Dependencies

Required packages:
```r
install.packages(c("Seurat", "Matrix", "glmnet", "caret", "igraph", 
                   "pROC", "ggplot2", "data.table", "dplyr", "future"))
```

Optional packages for additional features:
```r
install.packages(c("randomForest", "ranger", "e1071", "xgboost", 
                   "ComplexHeatmap", "ggalluvial", "patchwork"))
```

## 🔧 Quick Start

### Basic Assessment

```r
library(scClustEval)
library(Seurat)

# Load your Seurat object
seurat_obj <- readRDS("your_data.rds")

# Quick assessment
result <- RunAssessment(seurat_obj, cluster_col = "seurat_clusters")
print(result)

# Plot ROC curves
plot_roc(result)
```

### Clustering Optimization

```r
# Start with over-clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)

# Optimize clustering
seurat_obj <- RunOptimization(
  seurat_obj,
  cluster_col = "seurat_clusters",
  min_accuracy = 0.9,
  result_col = "optimized_clusters"
)

# Visualize results
DimPlot(seurat_obj, group.by = c("seurat_clusters", "optimized_clusters"))
```

### Using Expression Matrix Directly

```r
# Assessment with matrix
result <- sc_assessment(
  X = expression_matrix,  # cells x features
  labels = cluster_labels,
  classifier = "LR",
  n_per_class = 100
)

# Full optimization
optim_result <- sc_optimize_all(
  X = expression_matrix,
  labels = initial_clusters,
  min_accuracy = 0.9
)

final_clusters <- optim_result$final_labels
```

## 📊 Visualization

```r
# ROC curves
plot_roc(result, plot_type = "both")

# Confusion matrix heatmap
plot_confusion_heatmap(result, normalized = "R1")

# Optimization history
plot_optimization_history(optim_result)

# Cluster changes (Sankey diagram)
plot_cluster_sankey(initial_labels, final_labels)
```

## 🔬 Classifier Options

| Classifier | Code | Package |
|------------|------|---------|
| Logistic Regression | `"LR"` | glmnet |
| Random Forest | `"RF"` | randomForest |
| Fast Random Forest | `"RANGER"` | ranger |
| SVM | `"SVM"` | e1071 |
| Naive Bayes | `"NB"` | e1071 |
| Decision Tree | `"DT"` | rpart |
| XGBoost | `"XGB"` | xgboost |

## 📖 Algorithm

The core algorithm works by:

1. **Self-Projection**: Train a classifier to distinguish between clusters using a portion of the data
2. **Evaluation**: Assess prediction accuracy on held-out data
3. **Confusion Analysis**: Identify cluster pairs that are frequently confused
4. **Merging**: Merge clusters that cannot be well discriminated
5. **Iteration**: Repeat until target accuracy is reached

## 🤝 Citation

If you use scClustEval in your research, please cite:

```
Liu Z (2026). scClustEval: Single Cell Clustering Evaluation and Optimization 
Framework. R package version 1.0.0. https://github.com/Zaoqu-Liu/scClustEval
```

Also consider citing the original SCCAF paper:
```
Miao, Z., et al. (2020). Putative cell type discovery from single-cell gene
expression data. Nature Methods.
```

## 📬 Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: https://github.com/Zaoqu-Liu/scClustEval
- **Issues**: https://github.com/Zaoqu-Liu/scClustEval/issues

## 📄 License

MIT License © 2026 Zaoqu Liu
