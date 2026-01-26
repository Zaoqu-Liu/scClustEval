#' @title Visualization Functions for scClustEval
#' @description Functions for plotting assessment and optimization results
#' @name visualization
#'
#' @importFrom rlang .data
NULL

# Global variable declaration to avoid R CMD check notes
utils::globalVariables(c("stratum"))

#' @title Plot ROC curves
#' @description Plot ROC and Precision-Recall curves for assessment results
#'
#' @param result scClustEval result object from self_projection/sc_assessment
#' @param plot_type Type of plot: "both", "roc", or "prc"
#' @param show_auc Show AUC values on plot
#' @param show_cv Show cross-validation accuracy
#' @param show_acc Show test accuracy
#' @param colors Custom color palette
#' @param title Plot title
#' @param legend_position Legend position: "right", "bottom", or "none"
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' result <- sc_assessment(X, labels)
#' plot_roc(result)
#' plot_roc(result, plot_type = "roc")
#' }
#'
#' @import ggplot2
#' @export
plot_roc <- function(result,
                     plot_type = "both",
                     show_auc = TRUE,
                     show_cv = TRUE,
                     show_acc = TRUE,
                     colors = NULL,
                     title = NULL,
                     legend_position = "right") {

    if (!inherits(result, "scClustEval")) {
        stop("result must be an scClustEval object")
    }

    y_prob <- result$y_prob
    y_test <- result$y_test
    classes <- result$classes
    n_classes <- length(classes)

    # Set colors
    if (is.null(colors)) {
        colors <- .get_color_palette(n_classes)
    }

    # Calculate ROC and PRC for each class
    roc_data <- list()
    prc_data <- list()
    roc_aucs <- numeric(n_classes)
    prc_aucs <- numeric(n_classes)

    for (i in seq_along(classes)) {
        cl <- classes[i]

        # Binary labels for this class
        y_binary <- as.integer(y_test == cl)

        if (cl %in% colnames(y_prob)) {
            y_score <- y_prob[, cl]
        } else {
            warning(sprintf("Class %s not found in prediction probabilities", cl))
            next
        }

        # ROC curve using pROC
        roc_obj <- tryCatch({
            pROC::roc(y_binary, y_score, quiet = TRUE)
        }, error = function(e) NULL)

        if (!is.null(roc_obj)) {
            roc_data[[i]] <- data.frame(
                FPR = 1 - roc_obj$specificities,
                TPR = roc_obj$sensitivities,
                class = cl,
                stringsAsFactors = FALSE
            )
            roc_aucs[i] <- pROC::auc(roc_obj)

            # Precision-Recall
            # Sort by score
            ord <- order(y_score, decreasing = TRUE)
            y_sorted <- y_binary[ord]
            scores_sorted <- y_score[ord]

            # Calculate precision and recall at each threshold
            tp_cumsum <- cumsum(y_sorted)
            fp_cumsum <- cumsum(1 - y_sorted)
            total_positives <- sum(y_binary)

            # Handle edge case: no positive samples
            if (total_positives == 0) {
                warning(sprintf("No positive samples for class %s, skipping PRC", cl))
                next
            }

            precision <- tp_cumsum / (tp_cumsum + fp_cumsum)
            recall <- tp_cumsum / total_positives

            # Add endpoints
            precision <- c(1, precision)
            recall <- c(0, recall)

            prc_data[[i]] <- data.frame(
                Recall = recall,
                Precision = precision,
                class = cl,
                stringsAsFactors = FALSE
            )

            # AUC for PRC (trapezoidal)
            prc_aucs[i] <- tryCatch({
                sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
            }, error = function(e) NA)
        }
    }

    # Combine data
    roc_df <- do.call(rbind, roc_data)
    prc_df <- do.call(rbind, prc_data)

    names(roc_aucs) <- classes
    names(prc_aucs) <- classes

    # Min/Max AUC
    valid_roc_aucs <- roc_aucs[!is.na(roc_aucs)]
    valid_prc_aucs <- prc_aucs[!is.na(prc_aucs)]

    min_roc_auc <- min(valid_roc_aucs, na.rm = TRUE)
    max_roc_auc <- max(valid_roc_aucs, na.rm = TRUE)
    min_prc_auc <- min(valid_prc_aucs, na.rm = TRUE)
    max_prc_auc <- max(valid_prc_aucs, na.rm = TRUE)

    # Create annotation text
    make_annotation <- function(min_auc, max_auc, roc = TRUE) {
        text <- sprintf("AUC min: %.3f\nAUC max: %.3f", min_auc, max_auc)
        if (show_cv && !is.na(result$cv_accuracy)) {
            text <- paste0(text, sprintf("\nCV: %.3f", result$cv_accuracy))
        }
        if (show_acc) {
            text <- paste0(text, sprintf("\nTest: %.3f", result$accuracy))
        }
        text
    }

    # Build plots
    theme_roc <- theme_minimal() +
        theme(
            legend.position = legend_position,
            legend.title = element_blank(),
            panel.grid.minor = element_blank()
        )

    if (plot_type == "roc" || plot_type == "both") {
        p_roc <- ggplot(roc_df, aes(x = .data$FPR, y = .data$TPR, color = .data$class)) +
            geom_line(linewidth = 0.8) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
            scale_color_manual(values = colors) +
            labs(x = "False Positive Rate", y = "True Positive Rate") +
            coord_fixed() +
            xlim(0, 1) + ylim(0, 1) +
            theme_roc

        if (show_auc) {
            p_roc <- p_roc +
                annotate("text", x = 0.7, y = 0.2, hjust = 0,
                         label = make_annotation(min_roc_auc, max_roc_auc),
                         size = 3.5)
        }
    }

    if (plot_type == "prc" || plot_type == "both") {
        p_prc <- ggplot(prc_df, aes(x = .data$Recall, y = .data$Precision, color = .data$class)) +
            geom_line(linewidth = 0.8) +
            geom_abline(intercept = 1, slope = -1, linetype = "dashed", color = "gray50") +
            scale_color_manual(values = colors) +
            labs(x = "Recall", y = "Precision") +
            coord_fixed() +
            xlim(0, 1) + ylim(0, 1) +
            theme_roc

        if (show_auc) {
            p_prc <- p_prc +
                annotate("text", x = 0.05, y = 0.2, hjust = 0,
                         label = make_annotation(min_prc_auc, max_prc_auc, roc = FALSE),
                         size = 3.5)
        }
    }

    # Combine plots
    if (plot_type == "both") {
        if (requireNamespace("patchwork", quietly = TRUE)) {
            p <- p_roc + p_prc + patchwork::plot_layout(guides = "collect")
        } else {
            # Fallback without patchwork
            p <- p_roc
            message("Install 'patchwork' for side-by-side plots")
        }
    } else if (plot_type == "roc") {
        p <- p_roc
    } else {
        p <- p_prc
    }

    # Add title
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }

    return(p)
}

#' @title Plot confusion matrix heatmap
#' @description Visualize confusion matrix as a heatmap
#'
#' @param result scClustEval result object, or a confusion matrix directly
#' @param normalized Which matrix to plot: "raw", "R1", or "R2"
#' @param title Plot title
#' @param colors Color gradient (low to high)
#' @param show_values Show values in cells
#' @param text_size Size of text in cells
#'
#' @return A ggplot2 object
#'
#' @export
plot_confusion_heatmap <- function(result,
                                   normalized = "R1",
                                   title = NULL,
                                   colors = c("white", "gray20"),
                                   show_values = TRUE,
                                   text_size = 3) {

    # Extract matrix
    if (inherits(result, "scClustEval")) {
        if (normalized == "raw") {
            mat <- result$confusion_matrix
        } else if (normalized == "R1") {
            mat <- result$r1_normalized
        } else if (normalized == "R2") {
            mat <- result$r2_normalized
        } else {
            stop("normalized must be 'raw', 'R1', or 'R2'")
        }
    } else if (is.matrix(result)) {
        mat <- result
    } else {
        stop("result must be an scClustEval object or a matrix")
    }

    # Convert to long format
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("True", "Predicted", "Value")

    # Create plot
    p <- ggplot(df, aes(x = .data$Predicted, y = .data$True, fill = .data$Value)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = colors[1], high = colors[2]) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank()
        ) +
        labs(x = "Predicted", y = "True", fill = "Value")

    if (show_values) {
        # Format values
        df$label <- ifelse(df$Value > 0.01,
                           sprintf("%.2f", df$Value),
                           "")
        p <- p + geom_text(aes(label = .data$label),
                           size = text_size, color = "white")
    }

    if (!is.null(title)) {
        p <- p + ggtitle(title)
    } else {
        title_text <- switch(normalized,
                             "raw" = "Confusion Matrix",
                             "R1" = "R1-Normalized Confusion Matrix",
                             "R2" = "R2-Normalized Confusion Matrix")
        p <- p + ggtitle(title_text)
    }

    return(p)
}

#' @title Plot assessment summary
#' @description Create a comprehensive summary plot of assessment results
#'
#' @param result scClustEval result object
#' @param include Which plots to include: combination of "roc", "confusion", "accuracy"
#'
#' @return A combined ggplot2 object
#'
#' @export
plot_assessment_summary <- function(result,
                                    include = c("roc", "confusion", "accuracy")) {

    if (!inherits(result, "scClustEval")) {
        stop("result must be an scClustEval object")
    }

    plots <- list()

    if ("roc" %in% include) {
        plots$roc <- plot_roc(result, plot_type = "roc", legend_position = "none")
    }

    if ("confusion" %in% include) {
        plots$confusion <- plot_confusion_heatmap(result, normalized = "R1",
                                                  show_values = FALSE)
    }

    if ("accuracy" %in% include) {
        # Per-class accuracy bar plot
        acc_df <- data.frame(
            class = names(result$per_class_accuracy),
            accuracy = result$per_class_accuracy,
            stringsAsFactors = FALSE
        )
        acc_df <- acc_df[order(acc_df$accuracy), ]
        acc_df$class <- factor(acc_df$class, levels = acc_df$class)

        plots$accuracy <- ggplot(acc_df, aes(x = .data$class, y = .data$accuracy)) +
            geom_col(fill = "#3cb44b") +
            geom_hline(yintercept = result$accuracy, linetype = "dashed", color = "red") +
            coord_flip() +
            labs(x = NULL, y = "Accuracy", title = "Per-cluster Accuracy") +
            theme_minimal() +
            ylim(0, 1)
    }

    # Combine plots
    if (requireNamespace("patchwork", quietly = TRUE) && length(plots) > 1) {
        p <- patchwork::wrap_plots(plots, ncol = length(plots))
    } else {
        p <- plots[[1]]
    }

    return(p)
}

#' @title Plot optimization history
#' @description Visualize the progression of optimization rounds
#'
#' @param result scClustEval_optim result object from sc_optimize_all
#' @param metric Which metric to plot: "accuracy", "clusters", or "both"
#'
#' @return A ggplot2 object
#'
#' @export
plot_optimization_history <- function(result,
                                      metric = "both") {

    if (!inherits(result, "scClustEval_optim")) {
        stop("result must be an scClustEval_optim object")
    }

    # Create data frame
    n_rounds <- length(result$accuracy_history)

    df <- data.frame(
        Round = seq_len(n_rounds),
        Accuracy = result$accuracy_history,
        Clusters = result$n_clusters_history[-1],  # Remove initial
        stringsAsFactors = FALSE
    )

    if (metric == "accuracy") {
        p <- ggplot(df, aes(x = .data$Round, y = .data$Accuracy)) +
            geom_line(color = "#3cb44b", linewidth = 1) +
            geom_point(color = "#3cb44b", size = 3) +
            labs(x = "Round", y = "Accuracy", title = "Optimization Progress") +
            theme_minimal() +
            ylim(0, 1)

    } else if (metric == "clusters") {
        p <- ggplot(df, aes(x = .data$Round, y = .data$Clusters)) +
            geom_line(color = "#e6194b", linewidth = 1) +
            geom_point(color = "#e6194b", size = 3) +
            labs(x = "Round", y = "Number of Clusters", title = "Cluster Count") +
            theme_minimal()

    } else {
        # Both metrics
        df_long <- data.frame(
            Round = rep(df$Round, 2),
            Value = c(df$Accuracy, df$Clusters / max(df$Clusters)),
            Metric = rep(c("Accuracy", "Clusters (scaled)"), each = n_rounds),
            stringsAsFactors = FALSE
        )

        p <- ggplot(df_long, aes(x = .data$Round, y = .data$Value,
                                 color = .data$Metric, group = .data$Metric)) +
            geom_line(linewidth = 1) +
            geom_point(size = 3) +
            scale_color_manual(values = c("Accuracy" = "#3cb44b",
                                          "Clusters (scaled)" = "#e6194b")) +
            labs(x = "Round", y = "Value", title = "Optimization Progress",
                 color = NULL) +
            theme_minimal() +
            theme(legend.position = "bottom")

        # Add secondary axis
        p <- p + scale_y_continuous(
            name = "Accuracy",
            sec.axis = sec_axis(~ . * max(df$Clusters), name = "Clusters")
        )
    }

    return(p)
}

#' @title Plot Sankey diagram of cluster changes
#' @description Visualize cluster reassignment between rounds or annotations
#'
#' @param labels_from Original cluster labels
#' @param labels_to New cluster labels
#' @param title Plot title
#' @param colors Color palette
#' @param alpha Transparency of flows
#'
#' @return A ggplot2 object (using ggalluvial if available)
#'
#' @export
plot_cluster_sankey <- function(labels_from,
                                labels_to,
                                title = "Cluster Reassignment",
                                colors = NULL,
                                alpha = 0.6) {

    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
        stop("Package 'ggalluvial' is required for Sankey plots. ",
             "Install it with: install.packages('ggalluvial')")
    }

    # Create contingency table
    df <- data.frame(
        From = as.character(labels_from),
        To = as.character(labels_to),
        stringsAsFactors = FALSE
    )

    count_df <- as.data.frame(table(df$From, df$To))
    colnames(count_df) <- c("From", "To", "Freq")
    count_df <- count_df[count_df$Freq > 0, ]

    # Get all unique labels for colors
    all_labels <- unique(c(as.character(count_df$From), as.character(count_df$To)))

    if (is.null(colors)) {
        colors <- .get_color_palette(length(all_labels))
        names(colors) <- all_labels
    }

    # Create alluvial plot
    p <- ggplot(count_df,
                aes(y = .data$Freq,
                    axis1 = .data$From,
                    axis2 = .data$To)) +
        ggalluvial::geom_alluvium(aes(fill = .data$From), alpha = alpha) +
        ggalluvial::geom_stratum(width = 0.3, fill = "gray90", color = "gray50") +
        geom_text(stat = ggalluvial::StatStratum,
                  aes(label = after_stat(stratum)),
                  size = 3) +
        scale_x_discrete(limits = c("From", "To"), expand = c(0.15, 0.05)) +
        scale_fill_manual(values = colors) +
        labs(title = title, y = "Cell Count") +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            panel.grid = element_blank()
        )

    return(p)
}

#' @title Plot method for scClustEval objects
#' @param x scClustEval object
#' @param type Plot type: "roc", "confusion", "accuracy", "summary"
#' @param ... Additional arguments passed to specific plot functions
#' @export
#' @method plot scClustEval
plot.scClustEval <- function(x, type = "roc", ...) {
    if (type == "roc") {
        plot_roc(x, ...)
    } else if (type == "confusion") {
        plot_confusion_heatmap(x, ...)
    } else if (type == "summary") {
        plot_assessment_summary(x, ...)
    } else if (type == "accuracy") {
        plot_assessment_summary(x, include = "accuracy", ...)
    } else {
        stop("type must be 'roc', 'confusion', 'accuracy', or 'summary'")
    }
}

#' @title Plot cluster connections on embedding
#' @description Draw lines between cluster centroids on embedding space based on 
#'   confusion matrix values. This visualizes which clusters are frequently confused.
#'
#' @param embeddings Matrix of 2D embeddings (cells x 2), e.g., UMAP or t-SNE coordinates
#' @param labels Cluster labels for each cell
#' @param confusion_matrix Confusion/connection matrix between clusters (e.g., R1-normalized)
#' @param threshold Only draw lines for values above this threshold (default: 0)
#' @param line_color Color for connection lines (default: "#ffa500", orange)
#' @param line_scale Scale factor for line width (default: 10)
#' @param show_labels Show cluster labels at centroids (default: TRUE)
#' @param label_size Size of cluster labels (default: 4)
#' @param point_color Color for centroid points (default: "white")
#' @param point_size Size of centroid points (default: 5)
#' @param point_alpha Transparency of centroid points (default: 0.7)
#' @param add_to_plot If provided, add to existing ggplot object
#'
#' @return A ggplot2 object showing cluster connections
#'
#' @details
#' This function computes the median position (centroid) of each cluster in the
#' embedding space, then draws lines between centroids where the confusion matrix
#' value exceeds the threshold. Line width is proportional to the confusion value.
#'
#' @examples
#' \dontrun{
#' # Get UMAP coordinates from Seurat object
#' embeddings <- Seurat::Embeddings(seurat_obj, "umap")
#' labels <- seurat_obj$seurat_clusters
#'
#' # Run assessment and get R1 matrix
#' result <- sc_assessment(X, labels)
#'
#' # Plot connections
#' plot_cluster_links(embeddings, labels, result$r1_normalized, threshold = 0.1)
#' }
#'
#' @export
plot_cluster_links <- function(embeddings,
                               labels,
                               confusion_matrix,
                               threshold = 0,
                               line_color = "#ffa500",
                               line_scale = 10,
                               show_labels = TRUE,
                               label_size = 4,
                               point_color = "white",
                               point_size = 5,
                               point_alpha = 0.7,
                               add_to_plot = NULL) {
    
    # Validate inputs
    if (!is.matrix(embeddings) && !is.data.frame(embeddings)) {
        stop("embeddings must be a matrix or data.frame with 2 columns")
    }
    if (ncol(embeddings) < 2) {
        stop("embeddings must have at least 2 columns")
    }
    if (nrow(embeddings) != length(labels)) {
        stop("Number of rows in embeddings must match length of labels")
    }
    
    # Use first two columns
    embeddings <- as.data.frame(embeddings[, 1:2])
    colnames(embeddings) <- c("Dim1", "Dim2")
    labels <- as.factor(labels)
    
    # Compute centroids (median position for each cluster)
    centroids <- data.frame(
        cluster = levels(labels),
        Dim1 = numeric(nlevels(labels)),
        Dim2 = numeric(nlevels(labels))
    )
    
    for (i in seq_along(levels(labels))) {
        cl <- levels(labels)[i]
        mask <- labels == cl
        centroids$Dim1[i] <- median(embeddings$Dim1[mask])
        centroids$Dim2[i] <- median(embeddings$Dim2[mask])
    }
    rownames(centroids) <- centroids$cluster
    
    # Prepare confusion matrix
    cmat <- as.matrix(confusion_matrix)
    
    # Make symmetric and get upper triangle
    cmat_sym <- pmax(cmat, t(cmat), na.rm = TRUE)
    cmat_sym[lower.tri(cmat_sym, diag = TRUE)] <- 0
    
    # Find pairs above threshold
    pairs <- which(cmat_sym > threshold, arr.ind = TRUE)
    
    # Create line segments data
    if (nrow(pairs) > 0) {
        lines_df <- data.frame(
            x = numeric(nrow(pairs)),
            y = numeric(nrow(pairs)),
            xend = numeric(nrow(pairs)),
            yend = numeric(nrow(pairs)),
            weight = numeric(nrow(pairs))
        )
        
        cluster_names <- rownames(cmat)
        for (k in seq_len(nrow(pairs))) {
            i <- pairs[k, 1]
            j <- pairs[k, 2]
            
            cl_i <- cluster_names[i]
            cl_j <- cluster_names[j]
            
            if (cl_i %in% centroids$cluster && cl_j %in% centroids$cluster) {
                lines_df$x[k] <- centroids[cl_i, "Dim1"]
                lines_df$y[k] <- centroids[cl_i, "Dim2"]
                lines_df$xend[k] <- centroids[cl_j, "Dim1"]
                lines_df$yend[k] <- centroids[cl_j, "Dim2"]
                lines_df$weight[k] <- cmat_sym[i, j]
            }
        }
        
        # Remove rows with zero weight
        lines_df <- lines_df[lines_df$weight > 0, ]
    } else {
        lines_df <- data.frame(
            x = numeric(0), y = numeric(0),
            xend = numeric(0), yend = numeric(0),
            weight = numeric(0)
        )
    }
    
    # Create or add to plot
    if (is.null(add_to_plot)) {
        p <- ggplot()
    } else {
        p <- add_to_plot
    }
    
    # Add connection lines
    if (nrow(lines_df) > 0) {
        p <- p + geom_segment(
            data = lines_df,
            aes(x = .data$x, y = .data$y,
                xend = .data$xend, yend = .data$yend),
            linewidth = lines_df$weight * line_scale,
            color = line_color,
            alpha = 0.7,
            inherit.aes = FALSE
        )
    }
    
    # Add centroid points
    p <- p + geom_point(
        data = centroids,
        aes(x = .data$Dim1, y = .data$Dim2),
        color = point_color,
        size = point_size,
        alpha = point_alpha,
        inherit.aes = FALSE
    )
    
    # Add labels
    if (show_labels) {
        p <- p + geom_text(
            data = centroids,
            aes(x = .data$Dim1, y = .data$Dim2, label = .data$cluster),
            size = label_size,
            fontface = "bold",
            inherit.aes = FALSE
        )
    }
    
    # Theme
    if (is.null(add_to_plot)) {
        p <- p + theme_minimal() +
            labs(x = "Dimension 1", y = "Dimension 2") +
            coord_fixed()
    }
    
    return(p)
}

#' @title Plot cluster centroids on embedding
#' @description Mark cluster centroids on an embedding plot
#'
#' @param embeddings Matrix of 2D embeddings (cells x 2)
#' @param labels Cluster labels for each cell
#' @param point_color Color for centroid markers (default: "white")
#' @param point_size Size of centroid markers (default: 8)
#' @param point_alpha Transparency of markers (default: 0.6)
#' @param add_to_plot If provided, add to existing ggplot object
#'
#' @return A ggplot2 object with cluster centroids marked
#'
#' @export
plot_cluster_centers <- function(embeddings,
                                 labels,
                                 point_color = "white",
                                 point_size = 8,
                                 point_alpha = 0.6,
                                 add_to_plot = NULL) {
    
    # Use first two columns
    embeddings <- as.data.frame(embeddings[, 1:2])
    colnames(embeddings) <- c("Dim1", "Dim2")
    labels <- as.factor(labels)
    
    # Compute centroids
    centroids <- data.frame(
        cluster = levels(labels),
        Dim1 = numeric(nlevels(labels)),
        Dim2 = numeric(nlevels(labels))
    )
    
    for (i in seq_along(levels(labels))) {
        cl <- levels(labels)[i]
        mask <- labels == cl
        centroids$Dim1[i] <- median(embeddings$Dim1[mask])
        centroids$Dim2[i] <- median(embeddings$Dim2[mask])
    }
    
    # Create or add to plot
    if (is.null(add_to_plot)) {
        p <- ggplot()
    } else {
        p <- add_to_plot
    }
    
    # Add centroid markers
    p <- p + geom_point(
        data = centroids,
        aes(x = .data$Dim1, y = .data$Dim2),
        color = point_color,
        size = point_size,
        alpha = point_alpha,
        inherit.aes = FALSE
    )
    
    if (is.null(add_to_plot)) {
        p <- p + theme_minimal() +
            labs(x = "Dimension 1", y = "Dimension 2") +
            coord_fixed()
    }
    
    return(p)
}

#' @title Plot embedding with cluster links
#' @description Combined plot showing cells colored by cluster with confusion-based connections
#'
#' @param embeddings Matrix of 2D embeddings (cells x 2)
#' @param labels Cluster labels for each cell
#' @param confusion_matrix Confusion matrix between clusters
#' @param threshold Threshold for drawing connection lines (default: 0.1)
#' @param title Plot title
#' @param colors Optional color palette for clusters
#' @param point_size Size of cell points (default: 0.5)
#' @param line_color Color for connection lines (default: "#ffa500")
#' @param line_scale Scale factor for line widths (default: 10)
#' @param show_legend Show cluster legend (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @details
#' This function creates a scatter plot of cells colored by cluster assignment,
#' with lines connecting cluster centroids based on confusion matrix values.
#' Useful for visualizing which clusters are transcriptionally similar and
#' frequently confused by the classifier.
#'
#' @examples
#' \dontrun{
#' # With Seurat object
#' embeddings <- Seurat::Embeddings(seurat_obj, "umap")
#' labels <- seurat_obj$seurat_clusters
#' result <- sc_assessment(GetAssayData(seurat_obj), labels)
#'
#' plot_embedding_with_links(
#'   embeddings, labels, result$r1_normalized,
#'   threshold = 0.1, title = "UMAP with Confusion Links"
#' )
#' }
#'
#' @export
plot_embedding_with_links <- function(embeddings,
                                      labels,
                                      confusion_matrix,
                                      threshold = 0.1,
                                      title = NULL,
                                      colors = NULL,
                                      point_size = 0.5,
                                      line_color = "#ffa500",
                                      line_scale = 10,
                                      show_legend = TRUE) {
    
    # Prepare data
    embeddings <- as.data.frame(embeddings[, 1:2])
    colnames(embeddings) <- c("Dim1", "Dim2")
    labels <- as.factor(labels)
    
    # Set colors
    n_clusters <- nlevels(labels)
    if (is.null(colors)) {
        colors <- .get_color_palette(n_clusters)
    }
    
    # Create cell scatter data
    cell_df <- data.frame(
        Dim1 = embeddings$Dim1,
        Dim2 = embeddings$Dim2,
        Cluster = labels
    )
    
    # Base scatter plot
    p <- ggplot(cell_df, aes(x = .data$Dim1, y = .data$Dim2, color = .data$Cluster)) +
        geom_point(size = point_size, alpha = 0.6) +
        scale_color_manual(values = colors) +
        theme_minimal() +
        theme(
            legend.position = if (show_legend) "right" else "none"
        ) +
        labs(x = "Dimension 1", y = "Dimension 2", color = "Cluster")
    
    # Add connection lines
    p <- plot_cluster_links(
        embeddings = embeddings,
        labels = labels,
        confusion_matrix = confusion_matrix,
        threshold = threshold,
        line_color = line_color,
        line_scale = line_scale,
        show_labels = TRUE,
        label_size = 3,
        point_color = "white",
        point_size = 4,
        add_to_plot = p
    )
    
    # Add title
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    
    return(p)
}
