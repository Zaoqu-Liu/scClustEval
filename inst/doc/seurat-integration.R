## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  warning = FALSE,
  message = FALSE,
  eval = FALSE
)

## ----load, eval=FALSE---------------------------------------------------------
# library(scClustEval)
# library(ggplot2)

## ----load_seurat--------------------------------------------------------------
# library(Seurat)
# 
# # Load your preprocessed Seurat object
# # Assuming standard preprocessing: normalization, PCA, clustering
# seurat_obj <- readRDS("path/to/your/seurat_object.rds")
# 
# # Check available reductions and clusterings
# print(Reductions(seurat_obj))
# print(head(seurat_obj@meta.data))

## ----quick_assess-------------------------------------------------------------
# # Quick one-liner assessment
# accuracy <- QuickAssess(seurat_obj)
# 
# # Output:
# # Clustering Assessment
# # ---------------------
# # Test Accuracy:  85.3%
# # CV Accuracy:    84.7%
# # Max R1:         0.2341
# # Clusters:       15

## ----full_assess--------------------------------------------------------------
# # Run comprehensive assessment
# result <- RunAssessment(
#   object = seurat_obj,
#   cluster_col = "seurat_clusters",  # Which clustering to assess
#   use = "pca",                      # Use PCA embeddings
#   dims = 1:30,                      # First 30 PCs
#   classifier = "LR",                # Logistic regression
#   penalty = "l1",                   # L1 regularization
#   test_size = 0.5,
#   n_per_class = 100,
#   cv = 5,
#   seed = 42,
#   verbose = TRUE
# )
# 
# # View results
# print(result)
# 
# # Plot ROC curves
# plot_roc(result)

## ----feature_options----------------------------------------------------------
# # Option 1: PCA embeddings (recommended, default)
# result_pca <- RunAssessment(seurat_obj, use = "pca", dims = 1:30)
# 
# # Option 2: Normalized expression data
# result_data <- RunAssessment(seurat_obj, use = "data")
# 
# # Option 3: Highly variable genes only
# result_hvg <- RunAssessment(seurat_obj, use = "hvg")
# 
# # Option 4: Scaled data
# result_scaled <- RunAssessment(seurat_obj, use = "scale.data")
# 
# # Option 5: Other reductions (UMAP, Harmony, etc.)
# result_harmony <- RunAssessment(seurat_obj, use = "harmony", dims = 1:30)

## ----optimization-------------------------------------------------------------
# # Step 1: Create over-clustered data
# seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)
# initial_clusters <- length(unique(Idents(seurat_obj)))
# cat("Initial clusters:", initial_clusters, "\n")
# 
# # Step 2: Run optimization
# seurat_obj <- RunOptimization(
#   object = seurat_obj,
#   cluster_col = "seurat_clusters",     # Input clustering
#   result_col = "optimized_clusters",   # Output column name
#   prefix = "scClustEval",              # Prefix for intermediate columns
#   store_rounds = TRUE,                 # Store intermediate results
#   use = "pca",
#   dims = 1:30,
#   min_accuracy = 0.90,                 # Target 90% accuracy
#   max_rounds = 10,
#   classifier = "LR",
#   r1_cutoff = 0.5,
#   r2_cutoff = 0.05,
#   seed = 42,
#   verbose = TRUE
# )
# 
# # Step 3: Compare results
# final_clusters <- length(unique(seurat_obj$optimized_clusters))
# cat("Final clusters:", final_clusters, "\n")

## ----viz_optimization, fig.width=12-------------------------------------------
# # Compare clusterings side-by-side
# p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) +
#   labs(title = paste("Original (n =", initial_clusters, ")")) +
#   NoLegend()
# 
# p2 <- DimPlot(seurat_obj, group.by = "optimized_clusters", label = TRUE) +
#   labs(title = paste("Optimized (n =", final_clusters, ")")) +
#   NoLegend()
# 
# p1 + p2

## ----optim_history------------------------------------------------------------
# # Access stored optimization summary
# optim_summary <- seurat_obj@misc$scClustEval_optimization
# 
# # View metrics
# cat("Final accuracy:", optim_summary$final_accuracy, "\n")
# cat("Total rounds:", optim_summary$total_rounds, "\n")
# cat("Accuracy history:", optim_summary$accuracy_history, "\n")
# cat("Cluster count history:", optim_summary$n_clusters_history, "\n")

## ----constrained_optim--------------------------------------------------------
# # Create low and high resolution clusterings
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
# seurat_obj$low_res <- Idents(seurat_obj)
# 
# seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)
# seurat_obj$high_res <- Idents(seurat_obj)
# 
# # Optimize with constraint
# seurat_obj <- RunOptimization(
#   seurat_obj,
#   cluster_col = "high_res",
#   result_col = "constrained_clusters",
#   under_cluster_col = "low_res",  # Constraint: don't merge across these boundaries
#   min_accuracy = 0.95,
#   verbose = TRUE
# )
# 
# # Compare
# DimPlot(seurat_obj, group.by = c("low_res", "high_res", "constrained_clusters"),
#         ncol = 3, label = TRUE)

## ----reliability--------------------------------------------------------------
# # Run assessment first
# result <- RunAssessment(seurat_obj, cluster_col = "optimized_clusters")
# 
# # Add per-cluster reliability to metadata
# seurat_obj <- AddClusterReliability(
#   seurat_obj,
#   result = result,
#   cluster_col = "optimized_clusters",
#   reliability_col = "cluster_reliability"
# )
# 
# # Visualize reliability
# FeaturePlot(seurat_obj, features = "cluster_reliability",
#             cols = c("lightgray", "red"))
# 
# # Or as a violin plot per cluster
# VlnPlot(seurat_obj, features = "cluster_reliability",
#         group.by = "optimized_clusters", pt.size = 0)

## ----multi_assay--------------------------------------------------------------
# # Assess using RNA assay
# result_rna <- RunAssessment(
#   seurat_obj,
#   cluster_col = "seurat_clusters",
#   assay = "RNA",
#   use = "data"
# )
# 
# # Assess using integrated assay (if available)
# if ("integrated" %in% Assays(seurat_obj)) {
#   result_integrated <- RunAssessment(
#     seurat_obj,
#     cluster_col = "seurat_clusters",
#     assay = "integrated",
#     use = "pca"
#   )
# }

## ----v5_notes-----------------------------------------------------------------
# # scClustEval automatically handles V5 layer system
# # GetExpressionMatrix works for both V4 and V5
# 
# # Manual extraction if needed
# expr_matrix <- GetExpressionMatrix(
#   seurat_obj,
#   assay = "RNA",
#   slot = "data",       # "counts", "data", or "scale.data"
#   features = NULL      # NULL = all features
# )
# 
# # Check dimensions
# dim(expr_matrix)  # cells x features

## ----complete_workflow--------------------------------------------------------
# # ================================================
# # Complete scClustEval Workflow with Seurat
# # ================================================
# 
# # 1. Load and preprocess data
# seurat_obj <- readRDS("your_data.rds")
# 
# # 2. Initial assessment of default clustering
# cat("=== Initial Assessment ===\n")
# result_initial <- RunAssessment(
#   seurat_obj,
#   cluster_col = "seurat_clusters",
#   use = "pca",
#   dims = 1:30,
#   verbose = TRUE
# )
# cat("Initial accuracy:", round(result_initial$accuracy * 100, 1), "%\n")
# 
# # 3. Create over-clustered data
# seurat_obj <- FindClusters(seurat_obj, resolution = 2.0)
# seurat_obj$overclustered <- Idents(seurat_obj)
# 
# # 4. Optimize clustering
# cat("\n=== Optimization ===\n")
# seurat_obj <- RunOptimization(
#   seurat_obj,
#   cluster_col = "overclustered",
#   result_col = "final_clusters",
#   min_accuracy = 0.90,
#   max_rounds = 10,
#   verbose = TRUE
# )
# 
# # 5. Final assessment
# result_final <- RunAssessment(
#   seurat_obj,
#   cluster_col = "final_clusters",
#   use = "pca",
#   dims = 1:30,
#   verbose = FALSE
# )
# cat("\nFinal accuracy:", round(result_final$accuracy * 100, 1), "%\n")
# 
# # 6. Add reliability scores
# seurat_obj <- AddClusterReliability(seurat_obj, result_final,
#                                     cluster_col = "final_clusters")
# 
# # 7. Visualize
# p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) +
#   NoLegend() + ggtitle("Original")
# p2 <- DimPlot(seurat_obj, group.by = "final_clusters", label = TRUE) +
#   NoLegend() + ggtitle("Optimized")
# p3 <- FeaturePlot(seurat_obj, features = "cluster_reliability") +
#   ggtitle("Reliability")
# 
# p1 + p2 + p3
# 
# # 8. Save results
# saveRDS(seurat_obj, "seurat_with_optimized_clusters.rds")

## ----troubleshoot_pca---------------------------------------------------------
# # Solution: Run PCA first
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

## ----troubleshoot_small-------------------------------------------------------
# # glmnet may warn about small clusters
# # This is normal for clusters with < 8 cells
# # Consider filtering small clusters first
# cluster_sizes <- table(Idents(seurat_obj))
# keep_cells <- names(Idents(seurat_obj))[Idents(seurat_obj) %in% names(cluster_sizes[cluster_sizes >= 10])]
# seurat_subset <- subset(seurat_obj, cells = keep_cells)

## ----troubleshoot_memory------------------------------------------------------
# # Use PCA instead of full expression matrix
# result <- RunAssessment(seurat_obj, use = "pca", dims = 1:50)
# 
# # Or subsample
# result <- RunAssessment(seurat_obj, n_per_class = 50)  # Max 50 cells/cluster

