# scClustEval <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-Universe](https://zaoqu-liu.r-universe.dev/badges/scClustEval)](https://zaoqu-liu.r-universe.dev/scClustEval)
[![R-CMD-check](https://github.com/Zaoqu-Liu/scClustEval/workflows/R-CMD-check/badge.svg)](https://github.com/Zaoqu-Liu/scClustEval/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/scClustEval)](https://CRAN.R-project.org/package=scClustEval)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Single Cell Clustering Evaluation and Optimization Framework

**scClustEval** is a comprehensive R package designed for rigorous evaluation and iterative optimization of cell clustering in single-cell RNA sequencing (scRNA-seq) data. The package implements a self-projection machine learning framework that systematically assesses clustering reliability and identifies biologically indistinguishable cell populations for potential merging.

This package provides an R implementation inspired by the SCCAF (Single Cell Clustering Assessment Framework) Python package ([Miao et al., 2020, *Nature Methods*](https://doi.org/10.1038/s41592-020-0825-9)).

---

## Overview

Accurate cell type identification is fundamental to scRNA-seq analysis. However, conventional clustering algorithms often produce results that are sensitive to parameter choices and may not reflect true biological distinctions. **scClustEval** addresses this challenge through:

1. **Quantitative Assessment**: Objectively measures clustering quality using cross-validated classification accuracy
2. **Confusion Matrix Analysis**: Identifies cluster pairs with high misclassification rates indicating potential over-clustering
3. **Iterative Optimization**: Systematically merges indistinguishable clusters until a target discrimination accuracy is achieved

---

## Installation

### From R-Universe (Recommended)

```r
install.packages("scClustEval", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# Install remotes if not available
if (!require("remotes")) install.packages("remotes")

# Install scClustEval from GitHub
remotes::install_github("Zaoqu-Liu/scClustEval")
```

### System Requirements

- R (≥ 4.0.0)
- C++ compiler with C++11 support (for Rcpp components)

### Dependencies

**Core dependencies** (installed automatically):
- Seurat (≥ 4.0.0), SeuratObject, Matrix
- glmnet, caret, rpart, igraph
- pROC, ggplot2, rlang

**Optional dependencies** (for extended functionality):
```r
install.packages(c("randomForest", "ranger", "e1071", "xgboost", 
                   "leiden", "patchwork", "ggalluvial", "ComplexHeatmap"))
```

---

## Methodology

### Self-Projection Framework

The core algorithm employs a self-projection strategy:

1. **Data Partitioning**: Stratified splitting into training and test sets while preserving cluster proportions
2. **Classifier Training**: A multi-class classifier is trained to discriminate between clusters
3. **Cross-Validation**: Performance estimation via k-fold cross-validation on training data
4. **Confusion Matrix Computation**: Quantifies misclassification patterns between all cluster pairs
5. **Normalization**: Two complementary normalization schemes:
   - **R1**: Confusion rate relative to correctly classified cells
   - **R2**: Confusion rate relative to total cell count

### Optimization Strategy

The iterative optimization proceeds as follows:

1. Begin with over-clustered data (high resolution)
2. Assess clustering quality via self-projection
3. Identify cluster pairs exceeding confusion thresholds
4. Construct adjacency graph from confusion matrix
5. Apply community detection (Louvain/Leiden) to merge confused clusters
6. Repeat until target accuracy is achieved or convergence

---

## Usage

### Quick Assessment with Seurat Objects

```r
library(scClustEval)
library(Seurat)

# Load preprocessed Seurat object
seurat_obj <- readRDS("your_seurat_object.rds")

# Rapid assessment of current clustering
result <- RunAssessment(
  seurat_obj,
  cluster_col = "seurat_clusters",
  use = "pca",
  dims = 1:30,
  classifier = "LR"
)

# Summary statistics
print(result)
# Test Accuracy: 0.8532 (85.3%)
# CV Accuracy:   0.8467 (84.7%)
# Max R1:        0.2341
# Max R2:        0.0156

# Visualize ROC curves
plot_roc(result)
```

### Clustering Optimization Pipeline

```r
# Start with high-resolution clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)

# Iterative optimization
seurat_obj <- RunOptimization(
  seurat_obj,
  cluster_col = "seurat_clusters",
  result_col = "optimized_clusters",
  min_accuracy = 0.90,
  max_rounds = 10,
  classifier = "LR",
  r1_cutoff = 0.5,
  verbose = TRUE
)

# Compare clustering results
DimPlot(seurat_obj, group.by = c("seurat_clusters", "optimized_clusters"), 
        ncol = 2)
```

### Direct Matrix Input

For non-Seurat workflows:

```r
# Assessment with expression/embedding matrix
result <- sc_assessment(
  X = pca_embeddings,      # cells × features matrix
  labels = cluster_labels,
  classifier = "LR",
  penalty = "l1",
  test_size = 0.5,
  n_per_class = 100,
  cv = 5,
  seed = 42
)

# Full optimization pipeline
optim_result <- sc_optimize_all(
  X = pca_embeddings,
  labels = initial_clusters,
  min_accuracy = 0.90,
  r1_cutoff = 0.5,
  r2_cutoff = 0.05,
  classifier = "LR"
)

# Extract final clustering
final_clusters <- optim_result$final_labels
```

---

## Visualization Functions

```r
# ROC and Precision-Recall curves
plot_roc(result, plot_type = "both", show_auc = TRUE)

# R1-normalized confusion matrix heatmap
plot_confusion_heatmap(result, normalized = "R1")

# Optimization trajectory
plot_optimization_history(optim_result, metric = "both")

# Sankey diagram of cluster reassignments
plot_cluster_sankey(
  labels_from = initial_clusters,
  labels_to = final_clusters,
  title = "Cluster Optimization Flow"
)
```

---

## Supported Classifiers

| Classifier | Identifier | R Package | Notes |
|------------|------------|-----------|-------|
| Logistic Regression | `"LR"` | glmnet | L1/L2/Elastic-net regularization; **recommended** |
| Random Forest | `"RF"` | randomForest | Feature importance available |
| Ranger | `"RANGER"` | ranger | Fast RF implementation |
| Support Vector Machine | `"SVM"` | e1071 | RBF/linear/polynomial kernels |
| Naive Bayes | `"NB"` | e1071 | Efficient for high-dimensional data |
| Decision Tree | `"DT"` | rpart | Interpretable; feature importance |
| XGBoost | `"XGB"` | xgboost | Gradient boosting; requires installation |

---

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `classifier` | Machine learning algorithm | `"LR"` |
| `test_size` | Fraction of data for testing | `0.5` |
| `n_per_class` | Maximum training samples per cluster | `100` |
| `cv` | Cross-validation folds | `5` |
| `r1_cutoff` | R1 confusion threshold for merging | `0.5` |
| `r2_cutoff` | R2 confusion threshold for merging | `0.05` |
| `min_accuracy` | Target accuracy for optimization | `0.9` |
| `max_rounds` | Maximum optimization iterations | `10` |

---

## Performance

- **C++ Acceleration**: Core confusion matrix computations implemented in C++ via Rcpp/RcppArmadillo
- **Parallel Processing**: Optional multi-core support via the `future` framework
- **Memory Efficient**: Native sparse matrix support for large datasets
- **Seurat Compatible**: Full support for Seurat v4 and v5 object structures

---

## Citation

If you use **scClustEval** in your research, please cite:

```bibtex
@Manual{scClustEval2026,
  title = {scClustEval: Single Cell Clustering Evaluation and Optimization Framework},
  author = {Zaoqu Liu},
  year = {2026},
  note = {R package version 1.0.0},
  url = {https://github.com/Zaoqu-Liu/scClustEval}
}
```

Please also cite the original SCCAF methodology:

```bibtex
@Article{Miao2020,
  title = {Putative cell type discovery from single-cell gene expression data},
  author = {Miao, Zhichao and Moreno, Pablo and Huang, Ni and Papatheodorou, Irene and Brazma, Alvis and Teichmann, Sarah A.},
  journal = {Nature Methods},
  year = {2020},
  volume = {17},
  pages = {621--628},
  doi = {10.1038/s41592-020-0825-9}
}
```

---

## Related Resources

- **Original SCCAF**: [github.com/SCCAF/sccaf](https://github.com/SCCAF/sccaf) (Python implementation)
- **Publication**: [Miao et al., 2020, *Nature Methods*](https://doi.org/10.1038/s41592-020-0825-9)

---
## Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: [github.com/Zaoqu-Liu/scClustEval](https://github.com/Zaoqu-Liu/scClustEval)
- **Issues**: [github.com/Zaoqu-Liu/scClustEval/issues](https://github.com/Zaoqu-Liu/scClustEval/issues)

---

## License

MIT License © 2026 Zaoqu Liu

This package incorporates concepts from the SCCAF Python package, which is also released under the MIT License.
