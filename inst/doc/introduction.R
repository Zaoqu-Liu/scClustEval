## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
# # From GitHub
# devtools::install_github("Zaoqu-Liu/scClustEval")

## ----load, message=FALSE, warning=FALSE---------------------------------------
library(scClustEval)

## ----basic_assessment, eval=FALSE---------------------------------------------
# # Create example data
# set.seed(42)
# n_cells <- 500
# n_features <- 100
# n_clusters <- 5
# 
# # Generate expression matrix with cluster structure
# X <- matrix(0, nrow = n_cells, ncol = n_features)
# labels <- character(n_cells)
# 
# for (i in 1:n_clusters) {
#   idx <- ((i-1) * 100 + 1):(i * 100)
#   X[idx, ] <- matrix(rnorm(100 * n_features, mean = i), nrow = 100)
#   labels[idx] <- paste0("Cluster_", i)
# }
# 
# # Run assessment
# result <- sc_assessment(
#   X = X,
#   labels = labels,
#   classifier = "LR",
#   n_per_class = 50,
#   cv = 5
# )
# 
# # Print result
# print(result)

## ----seurat_example, eval=FALSE-----------------------------------------------
# library(Seurat)
# 
# # Load your Seurat object
# seurat_obj <- readRDS("your_data.rds")
# 
# # Run assessment on existing clustering
# result <- RunAssessment(
#   seurat_obj,
#   cluster_col = "seurat_clusters",
#   use = "pca",
#   dims = 1:30
# )
# 
# # View results
# print(result)
# 
# # Plot ROC curves
# plot_roc(result)

## ----optimization, eval=FALSE-------------------------------------------------
# # Start with over-clustering
# seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)
# 
# # Run optimization
# seurat_obj <- RunOptimization(
#   seurat_obj,
#   cluster_col = "seurat_clusters",
#   min_accuracy = 0.9,
#   result_col = "optimized_clusters"
# )
# 
# # Compare before and after
# DimPlot(seurat_obj, group.by = c("seurat_clusters", "optimized_clusters"))

## ----roc, eval=FALSE----------------------------------------------------------
# # Plot ROC and Precision-Recall curves
# plot_roc(result, plot_type = "both")
# 
# # ROC only
# plot_roc(result, plot_type = "roc")

## ----confusion, eval=FALSE----------------------------------------------------
# # Raw confusion matrix
# plot_confusion_heatmap(result, normalized = "raw")
# 
# # R1-normalized (used for merging decisions)
# plot_confusion_heatmap(result, normalized = "R1")

## ----history, eval=FALSE------------------------------------------------------
# # Run optimization with matrix input
# optim_result <- sc_optimize_all(
#   X = X,
#   labels = initial_labels,
#   min_accuracy = 0.9
# )
# 
# # Plot optimization progress
# plot_optimization_history(optim_result)

## ----classifiers, eval=FALSE--------------------------------------------------
# # Using different classifiers
# result_lr <- sc_assessment(X, labels, classifier = "LR")
# result_rf <- sc_assessment(X, labels, classifier = "RF")
# result_svm <- sc_assessment(X, labels, classifier = "SVM")

## ----constraints, eval=FALSE--------------------------------------------------
# # Create low and high resolution clusterings
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, key_added = "low_res")
# seurat_obj <- FindClusters(seurat_obj, resolution = 2.0, key_added = "high_res")
# 
# # Optimize with constraint
# seurat_obj <- RunOptimization(
#   seurat_obj,
#   cluster_col = "high_res",
#   under_cluster_col = "low_res",  # Constraint
#   min_accuracy = 0.95
# )

## ----parallel, eval=FALSE-----------------------------------------------------
# # Assessment uses parallel processing automatically
# # Control with n_cores parameter
# result <- sc_assessment(
#   X, labels,
#   n_cores = 4  # Use 4 cores
# )

## ----session------------------------------------------------------------------
sessionInfo()

