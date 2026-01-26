#' @title Utility Functions for scClustEval
#' @description Internal utility functions for data processing and manipulation
#' @name utils
NULL

#' @title Check Seurat availability
#' @description Internal function to check if Seurat package is available
#' @keywords internal
.check_seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for this function but not installed.\n",
             "Install with: install.packages('Seurat')\n",
             "Note: scClustEval does NOT auto-install Seurat to avoid version conflicts.",
             call. = FALSE)
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' is required for this function but not installed.\n",
             "Install with: install.packages('SeuratObject')",
             call. = FALSE)
    }
    TRUE
}

#' @title Make unique names
#' @description Make a name vector unique by adding suffix "_n"
#'
#' @param x Character vector with potential duplicates
#' @return Character vector with unique names
#'
#' @examples
#' make_unique_names(c("A", "B", "A", "C", "A"))
#' # Returns: c("A", "B", "A_1", "C", "A_2")
#'
#' @export
make_unique_names <- function(x) {
    x <- as.character(x)
    counts <- integer(length(x))
    seen <- new.env(hash = TRUE, parent = emptyenv())

    for (i in seq_along(x)) {
        name <- x[i]
        if (exists(name, envir = seen, inherits = FALSE)) {
            count <- get(name, envir = seen)
            x[i] <- paste0(name, "_", count)
            assign(name, count + 1L, envir = seen)
        } else {
            assign(name, 1L, envir = seen)
        }
    }
    return(x)
}

#' @title Stratified train-test split
#' @description Split data into training and test sets while maintaining class proportions
#'
#' @param X Feature matrix (cells x features)
#' @param y Class labels (factor or character vector)
#' @param test_size Fraction of data for testing (default: 0.5)
#' @param n_per_class Maximum number of samples per class in training set.
#'   If NULL, uses test_size fraction
#' @param seed Random seed for reproducibility
#'
#' @return A list with components:
#' \describe{
#'   \item{X_train}{Training feature matrix}
#'   \item{X_test}{Test feature matrix}
#'   \item{y_train}{Training labels}
#'   \item{y_test}{Test labels}
#'   \item{train_idx}{Indices of training samples}
#'   \item{test_idx}{Indices of test samples}
#' }
#'
#' @export
train_test_split_stratified <- function(X, y,
                                        test_size = 0.5,
                                        n_per_class = NULL,
                                        seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    y <- as.factor(y)
    classes <- levels(y)
    n_total <- length(y)

    train_idx <- integer(0)

    for (cl in classes) {
        cl_idx <- which(y == cl)
        n_cl <- length(cl_idx)

        if (!is.null(n_per_class)) {
            # Sample by number per class
            n_train <- min(n_per_class, floor(n_cl * (1 - test_size)))
            if (n_cl < 10) {
                n_train <- max(1, floor(n_cl * 0.9))
            }
        } else {
            # Sample by fraction
            n_train <- floor(n_cl * (1 - test_size))
        }

        n_train <- max(1, min(n_train, n_cl - 1))  # Ensure at least 1 in each set
        sampled <- sample(cl_idx, n_train)
        train_idx <- c(train_idx, sampled)
    }

    test_idx <- setdiff(seq_len(n_total), train_idx)

    # Subset data (works for both sparse and dense matrices)
    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]

    list(
        X_train = X_train,
        X_test = X_test,
        y_train = y[train_idx],
        y_test = y[test_idx],
        train_idx = train_idx,
        test_idx = test_idx
    )
}

#' @title Simple train-test split
#' @description Wrapper around train_test_split_stratified for simple usage
#'
#' @inheritParams train_test_split_stratified
#' @export
train_test_split <- train_test_split_stratified

#' @title Get expression matrix from Seurat object
#' @description Extract expression matrix from Seurat object with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param assay Assay to use (default: DefaultAssay)
#' @param slot Slot to use: "data", "counts", or "scale.data"
#' @param features Specific features to extract (default: all)
#'
#' @return Expression matrix (cells x features)
#'
#' @keywords internal
.get_expression_matrix <- function(object,
                                   assay = NULL,
                                   slot = "data",
                                   features = NULL) {
    if (!inherits(object, "Seurat")) {
        stop("Input must be a Seurat object")
    }

    if (is.null(assay)) {
        assay <- SeuratObject::DefaultAssay(object)
    }

    # Check Seurat version for compatibility
    seurat_version <- utils::packageVersion("Seurat")

    # Get the assay object
    assay_obj <- object[[assay]]

    # Extract data based on slot
    # Use GetAssayData which works for both V4 and V5
    mat <- tryCatch({
        SeuratObject::GetAssayData(object, slot = slot, assay = assay)
    }, error = function(e) {
        # Fallback for V5 layer-based access
        if (seurat_version >= "5.0.0") {
            layer_name <- switch(slot,
                                 "data" = "data",
                                 "counts" = "counts",
                                 "scale.data" = "scale.data",
                                 slot)
            # Try using Seurat's FetchData or direct access
            tryCatch({
                assay_obj[[layer_name]]
            }, error = function(e2) {
                stop(sprintf("Cannot extract %s slot from assay %s: %s",
                             slot, assay, e$message))
            })
        } else {
            stop(sprintf("Cannot extract %s slot: %s", slot, e$message))
        }
    })

    # Transpose to cells x features
    mat <- Matrix::t(mat)

    # Subset features if specified
    if (!is.null(features)) {
        features <- intersect(features, colnames(mat))
        if (length(features) == 0) {
            stop("None of the specified features found in the data")
        }
        mat <- mat[, features, drop = FALSE]
    }

    return(mat)
}

#' @title Get PCA embeddings from Seurat object
#' @description Extract PCA embeddings with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param reduction Name of reduction (default: "pca")
#' @param dims Dimensions to use (default: all)
#'
#' @return PCA embedding matrix (cells x dims)
#'
#' @keywords internal
.get_pca_embeddings <- function(object, reduction = "pca", dims = NULL) {
    if (!reduction %in% names(object@reductions)) {
        stop(sprintf("Reduction '%s' not found. Run PCA first.", reduction))
    }

    embeddings <- Seurat::Embeddings(object, reduction = reduction)

    if (!is.null(dims)) {
        if (max(dims) > ncol(embeddings)) {
            warning(sprintf("Requested %d dims but only %d available. Using all.",
                            max(dims), ncol(embeddings)))
            dims <- seq_len(ncol(embeddings))
        }
        embeddings <- embeddings[, dims, drop = FALSE]
    }

    return(embeddings)
}

#' @title Check and setup parallel backend
#' @description Configure parallel processing with future
#'
#' @param n_cores Number of cores to use. If NULL, uses all available
#' @param verbose Print information about parallel setup
#'
#' @return Number of cores configured
#'
#' @keywords internal
.setup_parallel <- function(n_cores = NULL, verbose = TRUE) {
    if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
    }
    n_cores <- max(1, as.integer(n_cores))

    # Use appropriate strategy based on OS
    os_type <- .Platform$OS.type

    if (n_cores == 1) {
        future::plan(future::sequential)
        if (verbose) message("Using sequential processing")
    } else {
        if (os_type == "unix") {
            # Unix/macOS: can use multicore (fork)
            future::plan(future::multicore, workers = n_cores)
        } else {
            # Windows: must use multisession
            future::plan(future::multisession, workers = n_cores)
        }
        if (verbose) message(sprintf("Using parallel processing with %d cores", n_cores))
    }

    return(n_cores)
}

#' @title Reset parallel backend
#' @description Reset to sequential processing
#' @keywords internal
.reset_parallel <- function() {
    future::plan(future::sequential)
}

#' @title Convert sparse matrix to dense if needed
#' @description Safely convert sparse matrix for operations that require dense
#'
#' @param mat Matrix (sparse or dense)
#' @param threshold Maximum number of elements before warning
#'
#' @return Dense matrix
#'
#' @keywords internal
.to_dense_matrix <- function(mat, threshold = 1e8) {
    if (!inherits(mat, "Matrix")) {
        return(as.matrix(mat))
    }

    n_elements <- nrow(mat) * ncol(mat)
    if (n_elements > threshold) {
        warning(sprintf(
            "Converting large sparse matrix (%d x %d = %.1e elements) to dense. This may use significant memory.",
            nrow(mat), ncol(mat), n_elements
        ))
    }

    return(as.matrix(mat))
}

#' @title Progress bar wrapper
#' @description Create a consistent progress bar
#'
#' @param n Total number of iterations
#' @param title Title for progress bar
#' @param verbose Whether to show progress
#'
#' @return Progress bar object or NULL
#'
#' @keywords internal
.create_progress_bar <- function(n, title = "", verbose = TRUE) {
    if (!verbose || !interactive()) return(NULL)

    pb <- txtProgressBar(min = 0, max = n, style = 3)
    return(pb)
}

#' @title Update progress bar
#' @keywords internal
.update_progress_bar <- function(pb, i) {
    if (!is.null(pb)) {
        setTxtProgressBar(pb, i)
    }
}

#' @title Close progress bar
#' @keywords internal
.close_progress_bar <- function(pb) {
    if (!is.null(pb)) {
        close(pb)
        cat("\n")
    }
}

#' @title Message with timestamp
#' @description Print message with timestamp for logging
#'
#' @param ... Message components
#' @param verbose Whether to print
#'
#' @keywords internal
.log_message <- function(..., verbose = TRUE) {
    if (verbose) {
        timestamp <- format(Sys.time(), "[%H:%M:%S]")
        message(timestamp, " ", ...)
    }
}

#' @title Validate input data
#' @description Check that input data is valid for assessment
#'
#' @param X Feature matrix
#' @param y Labels
#'
#' @keywords internal
.validate_input <- function(X, y) {
    # Check dimensions
    if (is.null(X) || is.null(y)) {
        stop("X and y cannot be NULL")
    }

    if (nrow(X) != length(y)) {
        stop(sprintf("Number of samples in X (%d) doesn't match length of y (%d)",
                     nrow(X), length(y)))
    }

    # Check for NA/NaN
    if (anyNA(y)) {
        stop("y contains NA values")
    }

    # Check number of classes
    n_classes <- length(unique(y))
    if (n_classes < 2) {
        stop("y must have at least 2 classes")
    }

    # Check minimum samples per class
    class_counts <- table(y)
    if (any(class_counts < 2)) {
        small_classes <- names(class_counts)[class_counts < 2]
        warning(sprintf("Classes with fewer than 2 samples: %s",
                        paste(small_classes, collapse = ", ")))
    }

    invisible(TRUE)
}

#' @title Get color palette
#' @description Get appropriate color palette for number of classes
#'
#' @param n Number of colors needed
#' @param palette Palette name: "default", "extended", or custom vector
#'
#' @return Character vector of colors
#'
#' @export
#' @keywords internal
.get_color_palette <- function(n, palette = "default") {
    if (is.character(palette) && length(palette) == 1) {
        if (palette == "default") {
            if (n <= 20) {
                colors <- .color_palette_20[seq_len(n)]
            } else {
                colors <- .color_palette_extended[seq_len(min(n, length(.color_palette_extended)))]
            }
        } else if (palette == "extended") {
            colors <- .color_palette_extended[seq_len(min(n, length(.color_palette_extended)))]
        } else {
            stop("Unknown palette name. Use 'default', 'extended', or provide a color vector.")
        }
    } else {
        colors <- palette[seq_len(min(n, length(palette)))]
    }

    # Extend with generated colors if needed
    if (length(colors) < n) {
        extra <- grDevices::colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c"))(n - length(colors))
        colors <- c(colors, extra)
    }

    return(colors)
}
