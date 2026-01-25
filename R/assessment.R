#' @title Assessment Functions for scClustEval
#' @description Core functions for clustering quality assessment using self-projection
#' @name assessment
NULL

#' @title Self-projection assessment
#' @description Core function for evaluating clustering quality using self-projection
#'
#' @param X Expression/feature matrix (cells x features). Can be sparse.
#' @param labels Cluster labels for each cell
#' @param classifier Classifier type: "LR", "RF", "SVM", "NB", "DT", "XGB", "RANGER"
#' @param penalty For LR: regularization type "l1", "l2", or "elasticnet"
#' @param lambda For LR: regularization strength. If NULL, uses CV to select
#' @param test_size Fraction of data for testing (default: 0.5)
#' @param n_per_class Maximum samples per class in training set. If NULL, uses test_size
#' @param cv Number of cross-validation folds on training set (0 to skip CV)
#' @param seed Random seed for reproducibility
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @param verbose Print progress messages
#'
#' @return A list (class "scClustEval") containing:
#' \describe{
#'   \item{accuracy}{Overall accuracy on test set}
#'   \item{cv_accuracy}{Mean cross-validation accuracy (if cv > 0)}
#'   \item{train_accuracy}{Accuracy on training set}
#'   \item{y_pred}{Predicted labels for test set}
#'   \item{y_test}{True labels for test set}
#'   \item{y_prob}{Prediction probabilities (matrix)}
#'   \item{confusion_matrix}{Confusion matrix}
#'   \item{r1_normalized}{R1-normalized confusion matrix}
#'   \item{r2_normalized}{R2-normalized confusion matrix}
#'   \item{per_class_accuracy}{Per-cluster accuracy}
#'   \item{classifier}{Trained classifier object}
#'   \item{classes}{Unique class labels}
#'   \item{max_r1}{Maximum R1 confusion value}
#'   \item{max_r2}{Maximum R2 confusion value}
#' }
#'
#' @details
#' The self-projection method works by:
#' 1. Splitting data into training and test sets (stratified by cluster)
#' 2. Training a classifier on the training set
#' 3. Evaluating prediction accuracy on the held-out test set
#' 4. Computing confusion matrices to identify poorly discriminated clusters
#'
#' High accuracy indicates that clusters are well-separated. Pairs of clusters
#' that are frequently confused may need to be merged.
#'
#' @examples
#' \dontrun{
#' # Basic usage with expression matrix
#' result <- self_projection(
#'   X = expression_matrix,
#'   labels = cluster_assignments,
#'   classifier = "LR"
#' )
#' print(result$accuracy)
#'
#' # With random forest
#' result <- self_projection(
#'   X = expression_matrix,
#'   labels = cluster_assignments,
#'   classifier = "RF",
#'   n_per_class = 100
#' )
#' }
#'
#' @seealso \code{\link{sc_assessment}}, \code{\link{plot_roc}}
#'
#' @export
self_projection <- function(X,
                            labels,
                            classifier = "LR",
                            penalty = "l1",
                            lambda = NULL,
                            test_size = 0.5,
                            n_per_class = NULL,
                            cv = 5,
                            seed = 1,
                            n_cores = NULL,
                            verbose = TRUE) {

    # Validate input
    .validate_input(X, labels)

    # Set seed
    if (!is.null(seed)) set.seed(seed)

    # Convert labels to factor
    labels <- as.factor(labels)
    classes <- sort(levels(labels))
    n_classes <- length(classes)

    if (verbose) {
        .log_message(sprintf("Starting self-projection with %d cells, %d features, %d classes",
                             nrow(X), ncol(X), n_classes))
    }

    # Split data
    split_data <- train_test_split_stratified(
        X = X,
        y = labels,
        test_size = test_size,
        n_per_class = n_per_class,
        seed = seed
    )

    X_train <- split_data$X_train
    X_test <- split_data$X_test
    y_train <- split_data$y_train
    y_test <- split_data$y_test

    if (verbose) {
        .log_message(sprintf("Split: %d training, %d test samples",
                             length(y_train), length(y_test)))
    }

    # Create and train classifier
    clf <- create_classifier(
        type = classifier,
        penalty = penalty,
        lambda = lambda,
        seed = seed
    )

    # Cross-validation on training set
    cv_accuracy <- NA
    if (cv > 0) {
        if (verbose) .log_message(sprintf("Running %d-fold cross-validation...", cv))
        cv_accuracy <- .run_cross_validation(
            X_train, y_train, classifier, penalty, lambda, cv, seed, n_cores
        )
        if (verbose) .log_message(sprintf("Mean CV accuracy: %.4f", cv_accuracy))
    }

    # Train on full training set
    if (verbose) .log_message("Training classifier...")
    clf$fit(X_train, y_train)

    # Evaluate on training set
    y_train_pred <- clf$predict(X_train)
    train_accuracy <- mean(y_train_pred == y_train)
    if (verbose) .log_message(sprintf("Training accuracy: %.4f", train_accuracy))

    # Evaluate on test set
    y_pred <- clf$predict(X_test)
    y_prob <- clf$predict_prob(X_test)
    test_accuracy <- mean(y_pred == y_test)
    if (verbose) .log_message(sprintf("Test accuracy: %.4f", test_accuracy))

    # Compute confusion matrix and normalizations
    cmat <- calc_confusion_matrix(y_test, y_pred, labels = classes)
    r1_norm <- normalize_confmat_r1(cmat)
    r2_norm <- normalize_confmat_r2(cmat)

    max_r1 <- max(r1_norm[upper.tri(r1_norm)], na.rm = TRUE)
    max_r2 <- max(r2_norm[upper.tri(r2_norm)], na.rm = TRUE)

    if (is.na(max_r1) || is.infinite(max_r1)) max_r1 <- 0
    if (is.na(max_r2) || is.infinite(max_r2)) max_r2 <- 0

    if (verbose) {
        .log_message(sprintf("Max R1: %.4f, Max R2: %.4f", max_r1, max_r2))
    }

    # Per-class accuracy
    class_acc <- per_cluster_accuracy(cmat)

    # Build result object
    result <- list(
        accuracy = test_accuracy,
        cv_accuracy = cv_accuracy,
        train_accuracy = train_accuracy,
        y_pred = y_pred,
        y_test = y_test,
        y_prob = y_prob,
        confusion_matrix = cmat,
        r1_normalized = r1_norm,
        r2_normalized = r2_norm,
        per_class_accuracy = class_acc,
        classifier = clf,
        classes = classes,
        max_r1 = max_r1,
        max_r2 = max_r2,
        n_train = length(y_train),
        n_test = length(y_test),
        call = match.call()
    )

    class(result) <- c("scClustEval", "list")
    return(result)
}

#' @title SCCAF Assessment
#' @description Alias for self_projection for compatibility
#'
#' @inheritParams self_projection
#' @inherit self_projection return
#'
#' @export
SCCAF_assessment <- function(...) {
    self_projection(...)
}

#' @title Single Cell Clustering Assessment
#' @description Main assessment function with user-friendly interface
#'
#' @inheritParams self_projection
#'
#' @return An scClustEval result object
#'
#' @examples
#' \dontrun{
#' # Assess clustering quality
#' result <- sc_assessment(
#'   X = expression_matrix,
#'   labels = seurat_object$seurat_clusters
#' )
#'
#' # Print summary
#' print(result)
#'
#' # Plot ROC curves
#' plot_roc(result)
#' }
#'
#' @seealso \code{\link{self_projection}}, \code{\link{RunAssessment}}
#'
#' @export
sc_assessment <- function(X,
                          labels,
                          classifier = "LR",
                          penalty = "l1",
                          lambda = NULL,
                          test_size = 0.5,
                          n_per_class = 100,
                          cv = 5,
                          seed = 1,
                          n_cores = NULL,
                          verbose = TRUE) {

    self_projection(
        X = X,
        labels = labels,
        classifier = classifier,
        penalty = penalty,
        lambda = lambda,
        test_size = test_size,
        n_per_class = n_per_class,
        cv = cv,
        seed = seed,
        n_cores = n_cores,
        verbose = verbose
    )
}

#' @title Run cross-validation
#' @description Internal function to perform stratified cross-validation
#'
#' @keywords internal
.run_cross_validation <- function(X, y, classifier, penalty, lambda, cv, seed, n_cores) {
    if (!is.null(seed)) set.seed(seed)

    y <- as.factor(y)
    n <- length(y)

    # Create stratified folds
    folds <- caret::createFolds(y, k = cv, returnTrain = FALSE)

    # Run CV in parallel or sequentially
    fold_accs <- numeric(cv)

    for (i in seq_len(cv)) {
        test_idx <- folds[[i]]
        train_idx <- setdiff(seq_len(n), test_idx)

        # Handle sparse matrices
        if (inherits(X, "Matrix")) {
            X_train_cv <- X[train_idx, , drop = FALSE]
            X_test_cv <- X[test_idx, , drop = FALSE]
        } else {
            X_train_cv <- X[train_idx, , drop = FALSE]
            X_test_cv <- X[test_idx, , drop = FALSE]
        }

        y_train_cv <- y[train_idx]
        y_test_cv <- y[test_idx]

        # Train and evaluate
        clf <- create_classifier(
            type = classifier,
            penalty = penalty,
            lambda = lambda,
            seed = seed
        )

        tryCatch({
            clf$fit(X_train_cv, y_train_cv)
            y_pred_cv <- clf$predict(X_test_cv)
            fold_accs[i] <- mean(y_pred_cv == y_test_cv)
        }, error = function(e) {
            warning(sprintf("CV fold %d failed: %s", i, e$message))
            fold_accs[i] <<- NA
        })
    }

    return(mean(fold_accs, na.rm = TRUE))
}

#' @title Get top markers from classifier
#' @description Extract top weighted features from a logistic regression classifier
#'
#' @param result scClustEval result object (from self_projection)
#' @param feature_names Names of features (genes). If NULL, uses column indices
#' @param top_n Number of top features per class (default: 10)
#'
#' @return Data frame with columns: class, feature, weight
#'
#' @export
get_top_markers <- function(result, feature_names = NULL, top_n = 10) {
    if (!inherits(result, "scClustEval")) {
        stop("result must be an scClustEval object")
    }

    clf <- result$classifier

    if (clf$type != "LR") {
        warning("get_top_markers is designed for logistic regression. ",
                "Using get_importance() for this classifier type.")

        imp <- clf$get_importance()
        if (is.null(imp)) {
            return(NULL)
        }

        # Format importance as data frame
        if (is.null(feature_names)) {
            feature_names <- names(imp)
            if (is.null(feature_names)) {
                feature_names <- seq_along(imp)
            }
        }

        df <- data.frame(
            class = "all",
            feature = feature_names,
            weight = as.numeric(imp),
            stringsAsFactors = FALSE
        )

        df <- df[order(-abs(df$weight)), ]
        df <- utils::head(df, top_n * length(result$classes))

        return(df)
    }

    # Get coefficients for LR
    coefs <- clf$get_coefficients()
    if (is.null(coefs)) {
        return(NULL)
    }

    # Process coefficients
    classes <- result$classes
    n_features <- length(coefs[[1]]) - 1  # Exclude intercept

    if (is.null(feature_names)) {
        feature_names <- paste0("Feature_", seq_len(n_features))
    }

    # Build data frame of weights
    all_weights <- list()

    for (i in seq_along(classes)) {
        cl <- classes[i]
        w <- as.numeric(coefs[[i]])[-1]  # Remove intercept

        # Get top positive weights for this class
        df <- data.frame(
            class = cl,
            feature = feature_names,
            weight = w,
            stringsAsFactors = FALSE
        )

        df <- df[df$weight > 0, ]
        df <- df[order(-df$weight), ]
        df <- utils::head(df, top_n)

        all_weights[[i]] <- df
    }

    result_df <- do.call(rbind, all_weights)
    rownames(result_df) <- NULL

    return(result_df)
}

#' @title Print method for scClustEval
#' @param x scClustEval object
#' @param ... Additional arguments (ignored)
#' @export
#' @method print scClustEval
print.scClustEval <- function(x, ...) {
    cat("\n============ scClustEval Assessment Result ============\n\n")

    cat(sprintf("  Test Accuracy:     %.4f (%.1f%%)\n",
                x$accuracy, x$accuracy * 100))

    if (!is.na(x$cv_accuracy)) {
        cat(sprintf("  CV Accuracy:       %.4f (%.1f%%)\n",
                    x$cv_accuracy, x$cv_accuracy * 100))
    }

    cat(sprintf("  Training Accuracy: %.4f (%.1f%%)\n",
                x$train_accuracy, x$train_accuracy * 100))

    cat(sprintf("\n  Samples: %d train, %d test\n", x$n_train, x$n_test))
    cat(sprintf("  Classes: %d (%s)\n",
                length(x$classes),
                paste(utils::head(x$classes, 5), collapse = ", ")))

    if (length(x$classes) > 5) {
        cat("           ...\n")
    }

    cat(sprintf("\n  Max R1 confusion: %.4f\n", x$max_r1))
    cat(sprintf("  Max R2 confusion: %.4f\n", x$max_r2))

    # Show lowest accuracy classes
    low_acc <- sort(x$per_class_accuracy)[1:min(3, length(x$per_class_accuracy))]
    cat("\n  Lowest per-class accuracies:\n")
    for (i in seq_along(low_acc)) {
        cat(sprintf("    %s: %.4f\n", names(low_acc)[i], low_acc[i]))
    }

    cat("\n=======================================================\n")

    invisible(x)
}

#' @title Summary method for scClustEval
#' @param object scClustEval object
#' @param ... Additional arguments (ignored)
#' @export
#' @method summary scClustEval
summary.scClustEval <- function(object, ...) {
    cat("\nscClustEval Assessment Summary\n")
    cat("------------------------------\n")

    # Basic metrics
    metrics <- data.frame(
        Metric = c("Test Accuracy", "CV Accuracy", "Train Accuracy",
                   "Max R1", "Max R2", "N Classes"),
        Value = c(
            sprintf("%.4f", object$accuracy),
            if (is.na(object$cv_accuracy)) "N/A" else sprintf("%.4f", object$cv_accuracy),
            sprintf("%.4f", object$train_accuracy),
            sprintf("%.4f", object$max_r1),
            sprintf("%.4f", object$max_r2),
            length(object$classes)
        ),
        stringsAsFactors = FALSE
    )
    print(metrics, row.names = FALSE)

    # Per-class summary
    cat("\nPer-class accuracy summary:\n")
    print(summary(object$per_class_accuracy))

    invisible(object)
}
