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
