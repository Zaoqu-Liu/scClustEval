#' @title Confusion Matrix Functions for scClustEval
#' @description Functions for computing and normalizing confusion matrices
#' @name confusion_matrix
NULL

#' @title Calculate confusion matrix
#' @description Compute confusion matrix from predictions and true labels
#'
#' @param y_true True class labels
#' @param y_pred Predicted class labels
#' @param labels Optional vector of labels to use (in order)
#'
#' @return A matrix with rows as true labels and columns as predicted labels
#'
#' @examples
#' y_true <- c("A", "A", "B", "B", "C", "C")
#' y_pred <- c("A", "B", "B", "B", "C", "A")
#' calc_confusion_matrix(y_true, y_pred)
#'
#' @export
calc_confusion_matrix <- function(y_true, y_pred, labels = NULL) {
    y_true <- as.factor(y_true)
    y_pred <- as.factor(y_pred)

    if (is.null(labels)) {
        labels <- sort(union(levels(y_true), levels(y_pred)))
    }

    # Ensure both factors have the same levels
    y_true <- factor(y_true, levels = labels)
    y_pred <- factor(y_pred, levels = labels)

    # Create confusion matrix
    cmat <- table(True = y_true, Predicted = y_pred)
    cmat <- as.matrix(cmat)

    return(cmat)
}

#' @title Normalize confusion matrix (R1 norm)
#' @description Normalize confusion matrix relative to correctly classified cells
#'
#' @details
#' R1 normalization measures the confusion rate between cluster pairs
#' relative to the number of correctly classified cells.
#'
#' For each pair (i, j), compute:
#' \deqn{R1(i,j) = \max\left(\frac{C_{ij}}{C_{jj}}, \frac{C_{ji}}{C_{ii}}\right)}
#'
#' where:
#' \itemize{
#'   \item \eqn{C_{ij}} = cells truly in cluster i but predicted as cluster j
#'   \item \eqn{C_{jj}} = cells truly in cluster j and correctly predicted (diagonal)
#' }
#'
#' The ratio \eqn{C_{ij}/C_{jj}} represents how many cells from cluster i
#' are misclassified as j, relative to the correctly classified cells in j.
#' A high R1 value indicates substantial confusion between the cluster pair.
#'
#' @param cmat Confusion matrix (from calc_confusion_matrix).
#'   Rows represent true labels, columns represent predicted labels.
#' @param mode Normalization mode: "1" (default, as in SCCAF) or "2"
#'
#' @return Symmetric matrix of pairwise R1-normalized confusion values.
#'   Values typically range from 0 to >1 (can exceed 1 when misclassifications
#'   outnumber correct classifications).
#'
#' @references
#' Miao, Z., et al. (2020). Putative cell type discovery from single-cell
#' gene expression data. Nature Methods.
#'
#' @examples
#' cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
#' rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")
#' normalize_confmat_r1(cmat)
#'
#' @export
normalize_confmat_r1 <- function(cmat, mode = "1") {
    cmat <- as.matrix(cmat)
    diag_vals <- diag(cmat)
    n <- nrow(cmat)

    # Initialize output matrix
    xmat <- matrix(0, nrow = n, ncol = n)
    rownames(xmat) <- rownames(cmat)
    colnames(xmat) <- colnames(cmat)

    for (i in seq_len(n)) {
        for (j in seq(i + 1, n)) {
            if (j > n) break

            if (mode == "1") {
                # Original SCCAF mode
                val1 <- if (diag_vals[j] > 0) cmat[i, j] / diag_vals[j] else 0
                val2 <- if (diag_vals[i] > 0) cmat[j, i] / diag_vals[i] else 0
            } else {
                # Alternative mode
                val1 <- if (diag_vals[i] > 0) cmat[i, j] / diag_vals[i] else 0
                val2 <- if (diag_vals[j] > 0) cmat[j, i] / diag_vals[j] else 0
            }

            xmat[i, j] <- xmat[j, i] <- max(val1, val2, na.rm = TRUE)
        }
    }

    return(xmat)
}

#' @title Normalize confusion matrix (R2 norm)
#' @description Normalize confusion matrix relative to total cell count
#'
#' @details
#' R2 normalization measures the overall impact of confusion between
#' cluster pairs on the entire dataset.
#'
#' For each pair (i, j), compute:
#' \deqn{R2(i,j) = \frac{C_{ij} + C_{ji}}{N}}
#'
#' where:
#' \itemize{
#'   \item \eqn{C_{ij}} = cells truly in cluster i but predicted as cluster j
#'   \item \eqn{C_{ji}} = cells truly in cluster j but predicted as cluster i
#'   \item \eqn{N} = total number of cells in the test set
#' }
#'
#' This gives the fraction of total cells that are confused between
#' clusters i and j. Unlike R1, R2 values are always between 0 and 1.
#'
#' @param cmat Confusion matrix (from calc_confusion_matrix).
#'   Rows represent true labels, columns represent predicted labels.
#'
#' @return Symmetric matrix of pairwise R2-normalized confusion values.
#'   Values range from 0 to 1, representing the fraction of total cells
#'   misclassified between each cluster pair.
#'
#' @references
#' Miao, Z., et al. (2020). Putative cell type discovery from single-cell
#' gene expression data. Nature Methods.
#'
#' @examples
#' cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
#' rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")
#' normalize_confmat_r2(cmat)
#'
#' @export
normalize_confmat_r2 <- function(cmat) {
    cmat <- as.matrix(cmat)
    n_total <- sum(cmat)

    if (n_total == 0) {
        warning("Total count is zero, returning zero matrix")
        return(cmat * 0)
    }

    # Sum symmetric off-diagonal elements
    xmat <- cmat + t(cmat)
    diag(xmat) <- 0

    # Normalize by total
    xmat <- xmat / n_total

    return(xmat)
}

#' @title Normalize confusion matrix
#' @description Wrapper function for confusion matrix normalization
#'
#' @param cmat Confusion matrix
#' @param method Normalization method: "R1" or "R2"
#' @param mode For R1: normalization mode "1" or "2"
#'
#' @return Normalized confusion matrix
#'
#' @export
normalize_confusion_matrix <- function(cmat, method = "R1", mode = "1") {
    method <- toupper(method)

    if (method == "R1") {
        return(normalize_confmat_r1(cmat, mode = mode))
    } else if (method == "R2") {
        return(normalize_confmat_r2(cmat))
    } else {
        stop("method must be 'R1' or 'R2'")
    }
}

#' @title Compute per-cluster accuracy
#' @description Calculate the classification accuracy for each cluster
#'
#' @param cmat Confusion matrix (rows = true, cols = predicted)
#'
#' @return Named numeric vector of per-cluster accuracies
#'
#' @examples
#' cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
#' rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")
#' per_cluster_accuracy(cmat)
#'
#' @export
per_cluster_accuracy <- function(cmat) {
    cmat <- as.matrix(cmat)
    diag_vals <- diag(cmat)
    row_sums <- rowSums(cmat)

    # Avoid division by zero
    accuracy <- ifelse(row_sums > 0, diag_vals / row_sums, 0)
    names(accuracy) <- rownames(cmat)

    return(accuracy)
}

#' @title Compute per-cell accuracy
#' @description Calculate the prediction confidence for each cell
#'
#' @param X Feature matrix (cells x features)
#' @param y True cell labels
#' @param clf Trained classifier object
#'
#' @return Numeric vector of per-cell accuracy scores
#'
#' @export
per_cell_accuracy <- function(X, y, clf) {
    y <- as.factor(y)

    # Get probability predictions
    probs <- clf$predict_prob(X)
    classes <- clf$get_classes()

    # Ensure column order matches
    if (!all(classes %in% colnames(probs))) {
        stop("Class names do not match between classifier and predictions")
    }

    # For each cell, get the probability of its true class
    n_cells <- nrow(X)
    accuracies <- numeric(n_cells)

    for (i in seq_len(n_cells)) {
        true_class <- as.character(y[i])
        if (true_class %in% colnames(probs)) {
            accuracies[i] <- probs[i, true_class]
        } else {
            accuracies[i] <- 0
        }
    }

    # Normalize by maximum probability
    max_probs <- apply(probs, 1, max)
    accuracies <- accuracies / pmax(max_probs, .Machine$double.eps)

    return(accuracies)
}

#' @title Get aggregated confusion statistics
#' @description Compute summary statistics from confusion matrix
#'
#' @param y_true True labels
#' @param y_pred Predicted labels
#' @param labels Optional ordered labels
#'
#' @return A list with confusion matrix and various statistics
#'
#' @keywords internal
.get_confusion_stats <- function(y_true, y_pred, labels = NULL) {
    cmat <- calc_confusion_matrix(y_true, y_pred, labels)
    r1mat <- normalize_confmat_r1(cmat)
    r2mat <- normalize_confmat_r2(cmat)

    # Overall accuracy
    accuracy <- sum(diag(cmat)) / sum(cmat)

    # Per-class stats
    per_class <- per_cluster_accuracy(cmat)

    # Maximum confusion values
    max_r1 <- max(r1mat[upper.tri(r1mat)], na.rm = TRUE)
    max_r2 <- max(r2mat[upper.tri(r2mat)], na.rm = TRUE)

    # Find most confused pairs
    r1_upper <- r1mat
    r1_upper[lower.tri(r1_upper, diag = TRUE)] <- 0
    max_idx <- which(r1_upper == max_r1, arr.ind = TRUE)

    confused_pairs <- if (nrow(max_idx) > 0) {
        data.frame(
            cluster1 = rownames(r1mat)[max_idx[, 1]],
            cluster2 = colnames(r1mat)[max_idx[, 2]],
            r1_value = max_r1,
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }

    list(
        confusion_matrix = cmat,
        r1_normalized = r1mat,
        r2_normalized = r2mat,
        accuracy = accuracy,
        per_class_accuracy = per_class,
        max_r1 = max_r1,
        max_r2 = max_r2,
        most_confused_pairs = confused_pairs
    )
}

#' @title Aggregate multiple confusion matrices
#' @description Take element-wise minimum across multiple confusion iterations
#'
#' @param cmat_list List of confusion matrices (from multiple sampling iterations)
#' @param method Aggregation method: "min" or "mean"
#'
#' @return Aggregated confusion matrix
#'
#' @keywords internal
.aggregate_confusion_matrices <- function(cmat_list, method = "min") {
    if (length(cmat_list) == 0) {
        stop("Empty confusion matrix list")
    }

    if (length(cmat_list) == 1) {
        return(cmat_list[[1]])
    }

    # Ensure all matrices have same dimensions
    dims <- sapply(cmat_list, dim)
    if (length(unique(dims[1, ])) > 1 || length(unique(dims[2, ])) > 1) {
        stop("All confusion matrices must have the same dimensions")
    }

    # Aggregate based on method
    if (method == "min") {
        result <- Reduce(pmin, cmat_list)
    } else if (method == "mean") {
        result <- Reduce(`+`, cmat_list) / length(cmat_list)
    } else {
        stop("method must be 'min' or 'mean'")
    }

    # Preserve names
    rownames(result) <- rownames(cmat_list[[1]])
    colnames(result) <- colnames(cmat_list[[1]])

    return(result)
}
