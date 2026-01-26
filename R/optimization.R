#' @title Optimization Functions for scClustEval
#' @description Functions for iterative clustering optimization
#' @name optimization
NULL

#' @title Single round of clustering optimization
#' @description Perform one round of confusion-matrix-guided cluster merging
#'
#' @param X Expression/feature matrix (cells x features)
#' @param labels Current cluster labels
#' @param classifier Classifier type: "LR", "RF", "SVM", etc.
#' @param penalty For LR: regularization type
#' @param lambda For LR: regularization strength
#' @param test_size Fraction for test set
#' @param n_per_class Max samples per class in training
#' @param cv Cross-validation folds (0 to skip)
#' @param n_iter Number of sampling iterations for confusion matrix (default: 3)
#' @param r1_cutoff Threshold for R1-normalized confusion (default: 0.1)
#' @param r2_cutoff Threshold for R2-normalized confusion (default: 0.05)
#' @param r1_mode R1 normalization mode: "1" or "2" (default: "1", as in SCCAF)
#' @param use_r1_only Use only R1 normalization for merging decisions
#' @param use_r2_only Use only R2 normalization for merging decisions
#' @param use_distance Use distance matrix in merging decision (default: FALSE)
#' @param dist_cutoff Distance cutoff for merging (default: 8.0)
#' @param use_projection Use self-projection labels for subsequent iterations (default: FALSE)
#' @param connection_matrix Optional connection matrix for constrained merging
#' @param resolution Louvain resolution for merging (default: 1.0)
#' @param seed Random seed
#' @param n_cores Number of cores for parallel processing
#' @param verbose Print progress
#'
#' @return A list containing:
#' \describe{
#'   \item{new_labels}{Merged cluster labels}
#'   \item{old_labels}{Original cluster labels}
#'   \item{group_map}{Mapping from old clusters to new groups}
#'   \item{accuracy}{Test accuracy from assessment}
#'   \item{cv_accuracy}{CV accuracy}
#'   \item{max_r1}{Maximum R1 confusion after this round}
#'   \item{max_r2}{Maximum R2 confusion after this round}
#'   \item{r1_matrix}{Aggregated R1-normalized confusion matrix}
#'   \item{r2_matrix}{Aggregated R2-normalized confusion matrix}
#'   \item{n_clusters_before}{Number of clusters before merging}
#'   \item{n_clusters_after}{Number of clusters after merging}
#'   \item{converged}{Whether optimization has converged (no merging possible)}
#'   \item{assessments}{List of assessment results from iterations}
#'   \item{self_projection_labels}{Self-projection predicted labels (if use_projection=TRUE)}
#' }
#'
#' @export
sc_optimize <- function(X,
                        labels,
                        classifier = "LR",
                        penalty = "l1",
                        lambda = NULL,
                        test_size = 0.5,
                        n_per_class = 100,
                        cv = 5,
                        n_iter = 3,
                        r1_cutoff = 0.1,
                        r2_cutoff = 0.05,
                        r1_mode = "1",
                        use_r1_only = FALSE,
                        use_r2_only = FALSE,
                        use_distance = FALSE,
                        dist_cutoff = 8.0,
                        use_projection = FALSE,
                        connection_matrix = NULL,
                        resolution = 1.0,
                        seed = 1,
                        n_cores = NULL,
                        verbose = TRUE) {

    # Convert labels
    labels <- as.factor(labels)
    old_labels <- labels
    classes <- sort(levels(labels))
    n_clusters_before <- length(classes)

    if (verbose) {
        .log_message(sprintf("Starting optimization round with %d clusters", n_clusters_before))
    }

    # Collect multiple iterations of confusion matrices
    r1_matrices <- list()
    r2_matrices <- list()
    cmatrices <- list()
    assessments <- list()
    accuracies <- numeric(n_iter)
    
    # Track labels for use_projection
    current_iter_labels <- labels
    self_proj_labels <- NULL

    for (i in seq_len(n_iter)) {
        iter_seed <- if (!is.null(seed)) seed + i - 1 else NULL

        if (verbose && i == 1) {
            .log_message("Running assessment iterations...")
        }

        # Run assessment - use projection labels for iterations after first if use_projection=TRUE
        train_labels <- if (i == 1 || !use_projection) labels else current_iter_labels
        
        result <- self_projection(
            X = X,
            labels = train_labels,
            classifier = classifier,
            penalty = penalty,
            lambda = lambda,
            test_size = test_size,
            n_per_class = n_per_class,
            cv = if (i == 1) cv else 0,  # Only CV on first iteration
            seed = iter_seed,
            n_cores = n_cores,
            verbose = (verbose && i == 1)
        )

        assessments[[i]] <- result
        accuracies[i] <- result$accuracy
        
        # Compute R1 with specified mode
        r1_matrices[[i]] <- normalize_confmat_r1(result$confusion_matrix, mode = r1_mode)
        r2_matrices[[i]] <- result$r2_normalized
        cmatrices[[i]] <- result$confusion_matrix
        
        # Store self-projection labels for use_projection
        if (i == 1) {
            # Get full predictions on all data for use_projection
            clf <- result$classifier
            all_pred <- clf$predict(X)
            self_proj_labels <- all_pred
            if (use_projection) {
                current_iter_labels <- as.factor(all_pred)
            }
        }
    }

    # Aggregate confusion matrices (take minimum across iterations)
    r1_mat <- .aggregate_confusion_matrices(r1_matrices, method = "min")
    r2_mat <- .aggregate_confusion_matrices(r2_matrices, method = "min")
    cmat_agg <- .aggregate_confusion_matrices(cmatrices, method = "min")

    max_r1 <- max(r1_mat[upper.tri(r1_mat)], na.rm = TRUE)
    max_r2 <- max(r2_mat[upper.tri(r2_mat)], na.rm = TRUE)

    if (is.na(max_r1) || is.infinite(max_r1)) max_r1 <- 0
    if (is.na(max_r2) || is.infinite(max_r2)) max_r2 <- 0

    if (verbose) {
        .log_message(sprintf("Aggregated Max R1: %.4f, Max R2: %.4f", max_r1, max_r2))
    }

    # Compute distance matrix if needed
    dist_mat <- NULL
    if (use_distance) {
        if (verbose) .log_message("Computing distance matrix...")
        dist_mat <- get_distance_matrix(X, labels, labels = classes, method = "euclidean")
        # Normalize: clusters are "close" if distance < cutoff
        dist_close <- dist_mat < dist_cutoff
    }

    # Determine which clusters to merge
    if (use_r1_only) {
        merge_mat <- r1_mat > r1_cutoff
    } else if (use_r2_only) {
        merge_mat <- r2_mat > r2_cutoff
    } else {
        # Use both R1 and R2 criteria (OR logic)
        if (!is.null(connection_matrix)) {
            # Constrained by connection matrix (AND logic)
            r1_exceed <- r1_mat > r1_cutoff
            merge_mat <- r1_exceed & (connection_matrix > 0)
        } else {
            merge_mat <- (r1_mat > r1_cutoff) | (r2_mat > r2_cutoff)
        }
        
        # Additionally constrain by distance if use_distance=TRUE
        if (use_distance) {
            # Only merge if confused AND close in distance
            merge_mat <- merge_mat & dist_close
        }
    }

    # Convert to numeric for clustering
    merge_mat <- merge_mat * 1L
    diag(merge_mat) <- 0

    # Cluster the merge matrix to get new groups
    groups <- cluster_adjacency_matrix(
        merge_mat,
        cutoff = 0,  # Already binarized
        resolution = resolution
    )

    # Check convergence
    n_new_clusters <- length(unique(groups))
    converged <- (n_new_clusters >= n_clusters_before)

    if (verbose) {
        if (converged) {
            .log_message("Converged: no clusters to merge")
        } else {
            .log_message(sprintf("Merging %d clusters into %d groups",
                                 n_clusters_before, n_new_clusters))
        }
    }

    # Create new labels - use projection labels if use_projection=TRUE
    source_labels <- if (use_projection && !is.null(self_proj_labels)) self_proj_labels else labels
    group_map <- setNames(as.character(groups), classes)
    new_labels <- group_map[as.character(source_labels)]
    new_labels <- factor(new_labels)

    # Results
    list(
        new_labels = new_labels,
        old_labels = old_labels,
        group_map = group_map,
        accuracy = min(accuracies),
        cv_accuracy = assessments[[1]]$cv_accuracy,
        max_r1 = max_r1,
        max_r2 = max_r2,
        r1_matrix = r1_mat,
        r2_matrix = r2_mat,
        confusion_matrix = cmat_agg,
        distance_matrix = dist_mat,
        n_clusters_before = n_clusters_before,
        n_clusters_after = n_new_clusters,
        converged = converged,
        assessments = assessments,
        self_projection_labels = self_proj_labels
    )
}

#' @title SCCAF Optimize
#' @description Alias for sc_optimize for compatibility
#' @inheritParams sc_optimize
#' @export
SCCAF_optimize <- sc_optimize

#' @title Full iterative optimization pipeline
#' @description Iteratively optimize clustering until target accuracy is reached
#'
#' @param X Expression/feature matrix (cells x features)
#' @param labels Initial cluster labels (should be over-clustered)
#' @param min_accuracy Target minimum accuracy to achieve (default: 0.9)
#' @param max_rounds Maximum optimization rounds (default: 10)
#' @param classifier Classifier type
#' @param penalty For LR: regularization type
#' @param lambda For LR: regularization strength
#' @param test_size Fraction for test set
#' @param n_per_class Max samples per class
#' @param cv CV folds
#' @param n_iter Iterations per round for confusion matrix
#' @param r1_cutoff Initial R1 cutoff (default: 0.5)
#' @param r2_cutoff Initial R2 cutoff (default: 0.05)
#' @param r1_step Step to reduce R1 cutoff each outer iteration (default: 0.01)
#' @param r2_step Step to reduce R2 cutoff each outer iteration (default: 0.001)
#' @param r1_mode R1 normalization mode: "1" or "2" (default: "1")
#' @param use_r1_only Use only R1 for merging
#' @param use_r2_only Use only R2 for merging
#' @param use_distance Use distance matrix in merging decisions (default: FALSE)
#' @param dist_cutoff Distance cutoff for merging (default: 8.0)
#' @param use_projection Use self-projection labels for subsequent iterations (default: FALSE)
#' @param under_cluster_labels Optional: under-clustering labels as constraint
#' @param min_outer_iter Minimum outer iterations before allowing convergence
#' @param seed Random seed
#' @param n_cores Number of cores
#' @param verbose Print progress
#'
#' @return A list containing:
#' \describe{
#'   \item{final_labels}{Final optimized cluster labels}
#'   \item{initial_labels}{Initial cluster labels}
#'   \item{round_history}{List of results from each round}
#'   \item{accuracy_history}{Vector of accuracies per round}
#'   \item{n_clusters_history}{Vector of cluster counts per round}
#'   \item{final_accuracy}{Final achieved accuracy}
#'   \item{total_rounds}{Number of rounds performed}
#'   \item{converged}{Whether target accuracy was reached}
#' }
#'
#' @details
#' The optimization proceeds in two levels:
#' \enumerate{
#'   \item \strong{Outer iterations}: Progressively lower the R1/R2 cutoffs
#'   \item \strong{Inner rounds}: Merge clusters based on current cutoffs
#' }
#'
#' The process continues until:
#' - Target accuracy is reached
#' - Maximum rounds are exceeded
#' - No more clusters can be merged
#'
#' @examples
#' \dontrun{
#' # Optimize over-clustered result
#' result <- sc_optimize_all(
#'   X = expression_matrix,
#'   labels = over_clustered_labels,
#'   min_accuracy = 0.9,
#'   classifier = "LR"
#' )
#'
#' # Get final labels
#' final_clusters <- result$final_labels
#' }
#'
#' @seealso \code{\link{sc_optimize}}, \code{\link{RunOptimization}}
#'
#' @export
sc_optimize_all <- function(X,
                            labels,
                            min_accuracy = 0.9,
                            max_rounds = 10,
                            classifier = "LR",
                            penalty = "l1",
                            lambda = NULL,
                            test_size = 0.5,
                            n_per_class = 100,
                            cv = 5,
                            n_iter = 3,
                            r1_cutoff = 0.5,
                            r2_cutoff = 0.05,
                            r1_step = 0.01,
                            r2_step = 0.001,
                            r1_mode = "1",
                            use_r1_only = FALSE,
                            use_r2_only = FALSE,
                            use_distance = FALSE,
                            dist_cutoff = 8.0,
                            use_projection = FALSE,
                            under_cluster_labels = NULL,
                            min_outer_iter = 3,
                            seed = 1,
                            n_cores = NULL,
                            verbose = TRUE) {

    # Initialize
    initial_labels <- as.factor(labels)
    current_labels <- initial_labels
    n_initial_clusters <- length(unique(current_labels))

    round_history <- list()
    accuracy_history <- numeric(0)
    n_clusters_history <- n_initial_clusters

    total_rounds <- 0
    current_accuracy <- 0
    outer_iter <- 0

    if (verbose) {
        .log_message(paste(rep("=", 60), collapse = ""))
        .log_message("Starting clustering optimization")
        .log_message(sprintf("Initial clusters: %d", n_initial_clusters))
        .log_message(sprintf("Target accuracy: %.2f%%", min_accuracy * 100))
        .log_message(paste(rep("=", 60), collapse = ""))
    }

    # Outer optimization loop
    while (current_accuracy < min_accuracy && outer_iter < 10) {
        outer_iter <- outer_iter + 1

        if (verbose) {
            .log_message(sprintf("\n--- Outer iteration %d ---", outer_iter))
            .log_message(sprintf("R1 cutoff: %.4f, R2 cutoff: %.4f", r1_cutoff, r2_cutoff))
        }

        # Inner optimization rounds
        for (round in seq_len(max_rounds)) {
            total_rounds <- total_rounds + 1

            if (verbose) {
                .log_message(sprintf("\nRound %d...", total_rounds))
            }

            # Prepare connection matrix if under-clustering constraint
            conn_mat <- NULL
            if (!is.null(under_cluster_labels)) {
                conn_mat <- get_connection_matrix(
                    under_cluster_labels,
                    current_labels
                )
            }

            # Run single optimization round
            round_result <- sc_optimize(
                X = X,
                labels = current_labels,
                classifier = classifier,
                penalty = penalty,
                lambda = lambda,
                test_size = test_size,
                n_per_class = n_per_class,
                cv = if (round == 1 && outer_iter == 1) cv else 0,
                n_iter = n_iter,
                r1_cutoff = r1_cutoff,
                r2_cutoff = r2_cutoff,
                r1_mode = r1_mode,
                use_r1_only = use_r1_only,
                use_r2_only = use_r2_only,
                use_distance = use_distance,
                dist_cutoff = dist_cutoff,
                use_projection = use_projection,
                connection_matrix = conn_mat,
                seed = if (!is.null(seed)) seed + total_rounds else NULL,
                n_cores = n_cores,
                verbose = verbose
            )

            # Store history
            round_history[[total_rounds]] <- round_result
            accuracy_history <- c(accuracy_history, round_result$accuracy)
            n_clusters_history <- c(n_clusters_history, round_result$n_clusters_after)

            current_accuracy <- round_result$accuracy

            if (verbose) {
                .log_message(sprintf("Accuracy: %.4f, Clusters: %d -> %d",
                                     current_accuracy,
                                     round_result$n_clusters_before,
                                     round_result$n_clusters_after))
            }

            # Check for convergence
            if (round_result$converged) {
                if (verbose) {
                    .log_message("Round converged (no clusters merged)")
                }
                break
            }

            # Check if target accuracy reached
            if (current_accuracy >= min_accuracy) {
                if (verbose) {
                    .log_message(sprintf("Target accuracy %.2f%% reached!",
                                         min_accuracy * 100))
                }
                break
            }

            # Check if only 1 cluster left
            if (round_result$n_clusters_after <= 1) {
                if (verbose) {
                    warning("Only 1 cluster remaining. Stopping optimization.")
                }
                break
            }

            # Update labels for next round
            current_labels <- round_result$new_labels
        }

        # Check overall convergence
        n_current <- length(unique(current_labels))
        if (current_accuracy >= min_accuracy) {
            break
        }

        # Check if cluster count is not decreasing across outer iterations
        # Compare with the last round's cluster count
        if (length(round_history) > 0) {
            last_n_clusters <- round_history[[length(round_history)]]$n_clusters_after
            if (n_current >= last_n_clusters && outer_iter >= min_outer_iter) {
                if (verbose) {
                    .log_message("Converged: cluster count not decreasing")
                }
                break
            }
        }

        # Update cutoffs for next outer iteration
        r1_cutoff <- round_result$max_r1 - r1_step
        r2_cutoff <- round_result$max_r2 - r2_step

        # Ensure cutoffs don't go negative
        r1_cutoff <- max(0.01, r1_cutoff)
        r2_cutoff <- max(0.001, r2_cutoff)
    }

    # Get final labels from last round
    if (length(round_history) > 0) {
        final_labels <- round_history[[length(round_history)]]$new_labels
    } else {
        final_labels <- current_labels
    }

    # Summary
    if (verbose) {
        .log_message(paste(rep("=", 60), collapse = ""))
        .log_message("Optimization complete")
        .log_message(sprintf("Final accuracy: %.4f (%.1f%%)",
                             current_accuracy, current_accuracy * 100))
        .log_message(sprintf("Clusters: %d -> %d",
                             n_initial_clusters, length(unique(final_labels))))
        .log_message(sprintf("Total rounds: %d", total_rounds))
        .log_message(paste(rep("=", 60), collapse = ""))
    }

    # Return results
    result <- list(
        final_labels = final_labels,
        initial_labels = initial_labels,
        round_history = round_history,
        accuracy_history = accuracy_history,
        n_clusters_history = n_clusters_history,
        final_accuracy = current_accuracy,
        total_rounds = total_rounds,
        converged = current_accuracy >= min_accuracy,
        call = match.call()
    )

    class(result) <- c("scClustEval_optim", "list")
    return(result)
}

#' @title SCCAF Optimize All
#' @description Alias for sc_optimize_all for compatibility
#' @inheritParams sc_optimize_all
#' @export
SCCAF_optimize_all <- sc_optimize_all

#' @title Print method for optimization result
#' @param x scClustEval_optim object
#' @param ... Additional arguments (ignored)
#' @export
#' @method print scClustEval_optim
print.scClustEval_optim <- function(x, ...) {
    cat("\n========== scClustEval Optimization Result ==========\n\n")

    cat(sprintf("  Final Accuracy: %.4f (%.1f%%)\n",
                x$final_accuracy, x$final_accuracy * 100))
    cat(sprintf("  Converged:      %s\n",
                if (x$converged) "Yes" else "No"))
    cat(sprintf("  Total Rounds:   %d\n", x$total_rounds))

    cat(sprintf("\n  Clusters: %d -> %d\n",
                length(unique(x$initial_labels)),
                length(unique(x$final_labels))))

    # Show accuracy progression
    if (length(x$accuracy_history) > 1) {
        cat("\n  Accuracy progression:\n")
        for (i in seq_along(x$accuracy_history)) {
            cat(sprintf("    Round %d: %.4f\n", i, x$accuracy_history[i]))
        }
    }

    cat("\n======================================================\n")

    invisible(x)
}

#' @title Summary method for optimization result
#' @param object scClustEval_optim object
#' @param ... Additional arguments (ignored)
#' @export
#' @method summary scClustEval_optim
summary.scClustEval_optim <- function(object, ...) {
    cat("\nscClustEval Optimization Summary\n")
    cat("--------------------------------\n")

    df <- data.frame(
        Round = seq_along(object$accuracy_history),
        Accuracy = object$accuracy_history,
        Clusters = object$n_clusters_history[-1],  # Remove initial
        stringsAsFactors = FALSE
    )

    print(df, row.names = FALSE)

    cat("\n")
    cat(sprintf("Initial clusters: %d\n", length(unique(object$initial_labels))))
    cat(sprintf("Final clusters:   %d\n", length(unique(object$final_labels))))
    cat(sprintf("Final accuracy:   %.4f\n", object$final_accuracy))

    invisible(object)
}
