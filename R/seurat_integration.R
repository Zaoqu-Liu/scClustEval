#' @title Seurat Integration Functions for scClustEval
#' @description Seamless integration with Seurat objects
#' @name seurat_integration
NULL

#' @title Run clustering assessment on Seurat object
#' @description Assess clustering quality directly on a Seurat object
#'
#' @param object Seurat object
#' @param cluster_col Column name in meta.data containing cluster labels
#'   (default: uses Idents)
#' @param assay Assay to use (default: DefaultAssay)
#' @param use Feature space to use: "raw" (normalized data), "pca", or
#'   other reduction name
#' @param dims Dimensions to use for PCA/reduction (default: 1:50)
#' @param features Optional: specific features to use. If NULL, uses all
#'   (for raw) or VariableFeatures
#' @param classifier Classifier type: "LR", "RF", etc.
#' @param penalty For LR: "l1" or "l2"
#' @param test_size Fraction for test set
#' @param n_per_class Max samples per class
#' @param cv Cross-validation folds
#' @param seed Random seed
#' @param n_cores Number of cores
#' @param verbose Print progress
#' @param ... Additional arguments passed to self_projection
#'
#' @return An scClustEval result object
#'
#' @details
#' This function extracts the appropriate data from the Seurat object
#' and runs self_projection assessment. By default, it uses PCA
#' coordinates if available, otherwise normalized data.
#'
#' For Seurat V4, uses GetAssayData with slots.
#' For Seurat V5, automatically uses LayerData when appropriate.
#'
#' @examples
#' \dontrun{
#' # Assess default clustering
#' result <- RunAssessment(seurat_obj)
#'
#' # Assess specific clustering with PCA
#' result <- RunAssessment(
#'   seurat_obj,
#'   cluster_col = "seurat_clusters",
#'   use = "pca",
#'   dims = 1:30
#' )
#'
#' # Assess with raw expression
#' result <- RunAssessment(
#'   seurat_obj,
#'   cluster_col = "manual_annotation",
#'   use = "raw"
#' )
#' }
#'
#' @seealso \code{\link{self_projection}}, \code{\link{RunOptimization}}
#'
#' @export
RunAssessment <- function(object,
                          cluster_col = NULL,
                          assay = NULL,
                          use = "pca",
                          dims = 1:50,
                          features = NULL,
                          classifier = "LR",
                          penalty = "l1",
                          test_size = 0.5,
                          n_per_class = 100,
                          cv = 5,
                          seed = 1,
                          n_cores = NULL,
                          verbose = TRUE,
                          ...) {

    # Validate input
    if (!inherits(object, "Seurat")) {
        stop("object must be a Seurat object")
    }

    # Get cluster labels
    if (is.null(cluster_col)) {
        labels <- Seurat::Idents(object)
    } else {
        if (!cluster_col %in% colnames(object@meta.data)) {
            stop(sprintf("Column '%s' not found in meta.data", cluster_col))
        }
        labels <- object@meta.data[[cluster_col]]
    }

    labels <- as.factor(labels)

    if (verbose) {
        .log_message(sprintf("Assessing clustering with %d cells, %d clusters",
                             length(labels), length(unique(labels))))
    }

    # Get feature matrix
    use <- tolower(use)

    if (use == "pca" || use %in% names(object@reductions)) {
        # Use reduction embeddings
        reduction_name <- if (use == "pca") "pca" else use

        if (!reduction_name %in% names(object@reductions)) {
            stop(sprintf("Reduction '%s' not found. Run Seurat::%s first.",
                         reduction_name,
                         if (reduction_name == "pca") "RunPCA" else "RunUMAP/etc"))
        }

        X <- .get_pca_embeddings(object, reduction = reduction_name, dims = dims)

        if (verbose) {
            .log_message(sprintf("Using %s embeddings (%d dimensions)",
                                 reduction_name, ncol(X)))
        }

    } else if (use == "raw" || use == "data") {
        # Use normalized expression data
        X <- .get_expression_matrix(object, assay = assay, slot = "data",
                                    features = features)

        if (verbose) {
            .log_message(sprintf("Using normalized data (%d features)", ncol(X)))
        }

    } else if (use == "counts") {
        X <- .get_expression_matrix(object, assay = assay, slot = "counts",
                                    features = features)

        if (verbose) {
            .log_message(sprintf("Using count data (%d features)", ncol(X)))
        }

    } else if (use == "scale.data" || use == "scaled") {
        X <- .get_expression_matrix(object, assay = assay, slot = "scale.data",
                                    features = features)

        if (verbose) {
            .log_message(sprintf("Using scaled data (%d features)", ncol(X)))
        }

    } else if (use == "hvg" || use == "variable") {
        # Use highly variable genes
        if (is.null(assay)) {
            assay <- SeuratObject::DefaultAssay(object)
        }

        var_features <- Seurat::VariableFeatures(object, assay = assay)
        if (length(var_features) == 0) {
            warning("No variable features found. Using all features.")
            var_features <- NULL
        }

        X <- .get_expression_matrix(object, assay = assay, slot = "data",
                                    features = var_features)

        if (verbose) {
            .log_message(sprintf("Using variable features (%d features)", ncol(X)))
        }

    } else {
        stop(sprintf("Unknown 'use' parameter: %s. Valid options: 'pca', 'raw', 'data', 'counts', 'scale.data', 'hvg'", use))
    }

    # Run assessment
    result <- self_projection(
        X = X,
        labels = labels,
        classifier = classifier,
        penalty = penalty,
        test_size = test_size,
        n_per_class = n_per_class,
        cv = cv,
        seed = seed,
        n_cores = n_cores,
        verbose = verbose,
        ...
    )

    # Add Seurat-specific info
    result$seurat_assay <- if (is.null(assay)) SeuratObject::DefaultAssay(object) else assay
    result$seurat_reduction <- if (use %in% c("pca", names(object@reductions))) use else NA
    result$seurat_cluster_col <- cluster_col

    return(result)
}

#' @title Run clustering optimization on Seurat object
#' @description Optimize clustering directly on a Seurat object
#'
#' @param object Seurat object
#' @param cluster_col Column name containing initial over-clustering
#' @param result_col Column name to store final optimized clustering
#' @param prefix Prefix for intermediate round columns (default: "scClustEval")
#' @param store_rounds Whether to store intermediate round results
#' @param assay Assay to use
#' @param use Feature space: "raw", "pca", etc.
#' @param dims Dimensions for PCA/reduction
#' @param features Specific features to use
#' @param min_accuracy Target accuracy
#' @param max_rounds Maximum rounds
#' @param classifier Classifier type
#' @param penalty For LR: regularization type
#' @param test_size Test fraction
#' @param n_per_class Max samples per class
#' @param cv CV folds
#' @param n_iter Confusion matrix iterations
#' @param r1_cutoff Initial R1 cutoff
#' @param r2_cutoff Initial R2 cutoff
#' @param under_cluster_col Optional: column with under-clustering constraint
#' @param seed Random seed
#' @param n_cores Number of cores
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return Modified Seurat object with optimized clustering added to meta.data
#'
#' @details
#' The function modifies the Seurat object by adding:
#' - Final optimized clustering in result_col (default: "scClustEval_clusters")
#' - Intermediate round results (if store_rounds = TRUE)
#' - Cluster reliability scores (optional)
#'
#' @examples
#' \dontrun{
#' # Basic optimization
#' seurat_obj <- RunOptimization(
#'   seurat_obj,
#'   cluster_col = "seurat_clusters",
#'   min_accuracy = 0.9
#' )
#'
#' # With under-clustering constraint
#' seurat_obj <- RunOptimization(
#'   seurat_obj,
#'   cluster_col = "high_res_clusters",
#'   under_cluster_col = "low_res_clusters",
#'   min_accuracy = 0.95
#' )
#' }
#'
#' @seealso \code{\link{sc_optimize_all}}, \code{\link{RunAssessment}}
#'
#' @export
RunOptimization <- function(object,
                            cluster_col,
                            result_col = "scClustEval_clusters",
                            prefix = "scClustEval",
                            store_rounds = TRUE,
                            assay = NULL,
                            use = "pca",
                            dims = 1:50,
                            features = NULL,
                            min_accuracy = 0.9,
                            max_rounds = 10,
                            classifier = "LR",
                            penalty = "l1",
                            test_size = 0.5,
                            n_per_class = 100,
                            cv = 5,
                            n_iter = 3,
                            r1_cutoff = 0.5,
                            r2_cutoff = 0.05,
                            under_cluster_col = NULL,
                            seed = 1,
                            n_cores = NULL,
                            verbose = TRUE,
                            ...) {

    # Validate input
    if (!inherits(object, "Seurat")) {
        stop("object must be a Seurat object")
    }

    if (!cluster_col %in% colnames(object@meta.data)) {
        stop(sprintf("Column '%s' not found in meta.data", cluster_col))
    }

    # Get initial labels
    labels <- as.factor(object@meta.data[[cluster_col]])

    # Get under-clustering labels if specified
    under_labels <- NULL
    if (!is.null(under_cluster_col)) {
        if (!under_cluster_col %in% colnames(object@meta.data)) {
            stop(sprintf("Under-clustering column '%s' not found", under_cluster_col))
        }
        under_labels <- as.factor(object@meta.data[[under_cluster_col]])
    }

    if (verbose) {
        .log_message(sprintf("Starting optimization from '%s' (%d clusters)",
                             cluster_col, length(unique(labels))))
    }

    # Get feature matrix (same logic as RunAssessment)
    use <- tolower(use)

    if (use == "pca" || use %in% names(object@reductions)) {
        reduction_name <- if (use == "pca") "pca" else use
        X <- .get_pca_embeddings(object, reduction = reduction_name, dims = dims)
    } else if (use == "raw" || use == "data") {
        X <- .get_expression_matrix(object, assay = assay, slot = "data",
                                    features = features)
    } else if (use == "hvg" || use == "variable") {
        if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
        var_features <- Seurat::VariableFeatures(object, assay = assay)
        if (length(var_features) == 0) var_features <- NULL
        X <- .get_expression_matrix(object, assay = assay, slot = "data",
                                    features = var_features)
    } else {
        X <- .get_expression_matrix(object, assay = assay, slot = "data")
    }

    # Store initial round
    if (store_rounds) {
        object@meta.data[[paste0(prefix, "_Round0")]] <- labels
    }

    # Run optimization
    optim_result <- sc_optimize_all(
        X = X,
        labels = labels,
        min_accuracy = min_accuracy,
        max_rounds = max_rounds,
        classifier = classifier,
        penalty = penalty,
        test_size = test_size,
        n_per_class = n_per_class,
        cv = cv,
        n_iter = n_iter,
        r1_cutoff = r1_cutoff,
        r2_cutoff = r2_cutoff,
        under_cluster_labels = under_labels,
        seed = seed,
        n_cores = n_cores,
        verbose = verbose,
        ...
    )

    # Store intermediate rounds
    if (store_rounds && length(optim_result$round_history) > 0) {
        for (i in seq_along(optim_result$round_history)) {
            round_labels <- optim_result$round_history[[i]]$new_labels
            col_name <- paste0(prefix, "_Round", i)
            object@meta.data[[col_name]] <- round_labels
        }
    }

    # Store final result
    object@meta.data[[result_col]] <- optim_result$final_labels

    # Store optimization summary as misc data
    object@misc[[paste0(prefix, "_optimization")]] <- list(
        final_accuracy = optim_result$final_accuracy,
        accuracy_history = optim_result$accuracy_history,
        n_clusters_history = optim_result$n_clusters_history,
        total_rounds = optim_result$total_rounds,
        converged = optim_result$converged
    )

    if (verbose) {
        .log_message(sprintf("Optimization complete. Final clusters stored in '%s'",
                             result_col))
    }

    return(object)
}

#' @title Add cluster reliability scores to Seurat object
#' @description Calculate and add per-cluster reliability scores
#'
#' @param object Seurat object
#' @param result scClustEval assessment result
#' @param cluster_col Cluster column that was assessed
#' @param reliability_col Name for reliability column
#'
#' @return Seurat object with reliability scores added
#'
#' @export
AddClusterReliability <- function(object,
                                  result,
                                  cluster_col = NULL,
                                  reliability_col = NULL) {

    if (!inherits(object, "Seurat")) {
        stop("object must be a Seurat object")
    }

    if (!inherits(result, "scClustEval")) {
        stop("result must be an scClustEval object")
    }

    # Get cluster labels
    if (is.null(cluster_col)) {
        cluster_col <- result$seurat_cluster_col
        if (is.null(cluster_col)) {
            cluster_col <- "seurat_clusters"
        }
    }

    if (!cluster_col %in% colnames(object@meta.data)) {
        stop(sprintf("Column '%s' not found in meta.data", cluster_col))
    }

    # Generate reliability column name
    if (is.null(reliability_col)) {
        reliability_col <- paste0(cluster_col, "_reliability")
    }

    # Get per-cluster accuracies
    cluster_acc <- result$per_class_accuracy

    # Map to cells
    labels <- as.character(object@meta.data[[cluster_col]])
    reliability <- cluster_acc[labels]

    # Handle any missing
    reliability[is.na(reliability)] <- 0

    # Add to object
    object@meta.data[[reliability_col]] <- reliability

    message(sprintf("Added reliability scores to '%s'", reliability_col))

    return(object)
}

#' @title Extract feature matrix from Seurat (exported helper)
#' @description Helper function to extract feature matrix with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param assay Assay name
#' @param slot Slot name: "data", "counts", "scale.data"
#' @param features Features to extract
#'
#' @return Matrix (cells x features)
#'
#' @export
GetExpressionMatrix <- function(object, assay = NULL, slot = "data", features = NULL) {
    .get_expression_matrix(object, assay = assay, slot = slot, features = features)
}

#' @title Quick assessment from Seurat object
#' @description One-liner for quick clustering quality check
#'
#' @param object Seurat object
#' @param ... Additional arguments passed to RunAssessment
#'
#' @return Prints summary and returns accuracy
#'
#' @examples
#' \dontrun{
#' # Quick check
#' accuracy <- QuickAssess(seurat_obj)
#' }
#'
#' @export
QuickAssess <- function(object, ...) {
    result <- RunAssessment(object, verbose = FALSE, ...)
    cat(sprintf("\nClustering Assessment\n"))
    cat(sprintf("---------------------\n"))
    cat(sprintf("Test Accuracy:  %.1f%%\n", result$accuracy * 100))
    if (!is.na(result$cv_accuracy)) {
        cat(sprintf("CV Accuracy:    %.1f%%\n", result$cv_accuracy * 100))
    }
    cat(sprintf("Max R1:         %.4f\n", result$max_r1))
    cat(sprintf("Clusters:       %d\n", length(result$classes)))
    invisible(result$accuracy)
}
