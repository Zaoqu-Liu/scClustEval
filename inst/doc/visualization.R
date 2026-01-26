## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  warning = FALSE,
  message = FALSE,
  dpi = 100,
  eval = FALSE
)
set.seed(42)

## ----load---------------------------------------------------------------------
# library(scClustEval)
# library(ggplot2)

## ----prepare_data-------------------------------------------------------------
# # Create synthetic data
# set.seed(42)
# n_cells <- 800
# n_features <- 60
# n_clusters <- 5
# 
# X <- matrix(nrow = n_cells, ncol = n_features)
# labels <- character(n_cells)
# cells_per_cluster <- n_cells / n_clusters
# 
# for (i in 1:n_clusters) {
#   start_idx <- (i - 1) * cells_per_cluster + 1
#   end_idx <- i * cells_per_cluster
#   cluster_mean <- rnorm(n_features, mean = i * 1.5, sd = 0.3)
#   X[start_idx:end_idx, ] <- matrix(
#     rep(cluster_mean, cells_per_cluster) + rnorm(cells_per_cluster * n_features, sd = 0.8),
#     nrow = cells_per_cluster,
#     byrow = TRUE
#   )
#   labels[start_idx:end_idx] <- paste0("Type_", LETTERS[i])
# }
# 
# colnames(X) <- paste0("Gene_", 1:n_features)
# 
# # Run assessment
# result <- sc_assessment(
#   X = X, labels = labels,
#   classifier = "LR",
#   n_per_class = 100,
#   cv = 5,
#   seed = 42,
#   verbose = FALSE
# )

## ----roc_basic----------------------------------------------------------------
# plot_roc(result, plot_type = "roc")

## ----prc_basic----------------------------------------------------------------
# plot_roc(result, plot_type = "prc")

## ----roc_both, fig.width=12, fig.height=5-------------------------------------
# plot_roc(result, plot_type = "both", show_auc = TRUE, show_cv = TRUE, show_acc = TRUE)

## ----roc_custom, fig.width=8--------------------------------------------------
# # Custom colors
# custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
# 
# plot_roc(
#   result,
#   plot_type = "roc",
#   colors = custom_colors,
#   title = "Cluster Discrimination Performance",
#   legend_position = "bottom"
# )

## ----confusion_raw------------------------------------------------------------
# plot_confusion_heatmap(result, normalized = "raw", title = "Raw Confusion Matrix")

## ----confusion_r1-------------------------------------------------------------
# plot_confusion_heatmap(
#   result,
#   normalized = "R1",
#   title = "R1-Normalized Confusion",
#   show_values = TRUE,
#   text_size = 4
# )

## ----confusion_r2-------------------------------------------------------------
# plot_confusion_heatmap(
#   result,
#   normalized = "R2",
#   title = "R2-Normalized Confusion"
# )

## ----confusion_colors---------------------------------------------------------
# # Custom gradient
# plot_confusion_heatmap(
#   result,
#   normalized = "R1",
#   colors = c("#F7FBFF", "#08306B"),  # Blue gradient
#   title = "Blue Theme Confusion Matrix"
# )

## ----confusion_compare, fig.width=14, fig.height=5----------------------------
# library(gridExtra)
# 
# p1 <- plot_confusion_heatmap(result, normalized = "raw", title = "Raw Counts")
# p2 <- plot_confusion_heatmap(result, normalized = "R1", title = "R1 Normalized")
# p3 <- plot_confusion_heatmap(result, normalized = "R2", title = "R2 Normalized")
# 
# grid.arrange(p1, p2, p3, ncol = 3)

## ----assessment_summary, fig.width=10-----------------------------------------
# plot_assessment_summary(result, include = c("accuracy"))

## ----accuracy_custom----------------------------------------------------------
# # Extract per-cluster accuracy
# acc_df <- data.frame(
#   Cluster = names(result$per_class_accuracy),
#   Accuracy = result$per_class_accuracy
# )
# acc_df <- acc_df[order(acc_df$Accuracy), ]
# acc_df$Cluster <- factor(acc_df$Cluster, levels = acc_df$Cluster)
# 
# ggplot(acc_df, aes(x = Cluster, y = Accuracy, fill = Accuracy)) +
#   geom_col(width = 0.7) +
#   geom_hline(yintercept = result$accuracy, linetype = "dashed", color = "red", size = 1) +
#   geom_text(aes(label = sprintf("%.1f%%", Accuracy * 100)),
#             hjust = -0.1, size = 3.5) +
#   scale_fill_gradient2(low = "#d62728", mid = "#ff7f0e", high = "#2ca02c",
#                        midpoint = 0.85, limits = c(0.7, 1)) +
#   coord_flip() +
#   labs(title = "Per-Cluster Classification Accuracy",
#        subtitle = sprintf("Overall accuracy: %.1f%% (dashed line)", result$accuracy * 100),
#        x = NULL, y = "Accuracy") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"),
#         legend.position = "none") +
#   ylim(0, 1.15)

## ----optim_data---------------------------------------------------------------
# # Create over-clustered scenario
# labels_over <- labels
# labels_over[labels == "Type_A"][1:80] <- "Type_A1"
# labels_over[labels == "Type_A"][81:160] <- "Type_A2"
# labels_over[labels == "Type_B"][1:80] <- "Type_B1"
# labels_over[labels == "Type_B"][81:160] <- "Type_B2"
# 
# # Run optimization
# optim_result <- sc_optimize_all(
#   X = X,
#   labels = labels_over,
#   min_accuracy = 0.90,
#   max_rounds = 8,
#   classifier = "LR",
#   r1_cutoff = 0.5,
#   seed = 42,
#   verbose = FALSE
# )

## ----optim_history_acc--------------------------------------------------------
# plot_optimization_history(optim_result, metric = "accuracy")

## ----optim_history_cluster----------------------------------------------------
# plot_optimization_history(optim_result, metric = "clusters")

## ----optim_history_both, fig.width=10-----------------------------------------
# plot_optimization_history(optim_result, metric = "both")

## ----optim_custom, fig.width=10-----------------------------------------------
# # Create detailed optimization trajectory
# rounds <- seq_along(optim_result$accuracy_history)
# 
# df_optim <- data.frame(
#   Round = rounds,
#   Accuracy = optim_result$accuracy_history,
#   Clusters = optim_result$n_clusters_history[-1]
# )
# 
# p1 <- ggplot(df_optim, aes(x = Round, y = Accuracy)) +
#   geom_ribbon(aes(ymin = 0.7, ymax = Accuracy), fill = "#3cb44b", alpha = 0.3) +
#   geom_line(color = "#3cb44b", size = 1.5) +
#   geom_point(color = "#3cb44b", size = 4) +
#   geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", size = 1) +
#   annotate("text", x = max(rounds) - 0.5, y = 0.92,
#            label = "Target", color = "red", fontface = "bold") +
#   scale_y_continuous(labels = scales::percent, limits = c(0.7, 1)) +
#   labs(title = "Accuracy Improvement", y = "Accuracy", x = "Round") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold", size = 14))
# 
# p2 <- ggplot(df_optim, aes(x = Round, y = Clusters)) +
#   geom_area(fill = "#e6194b", alpha = 0.3) +
#   geom_line(color = "#e6194b", size = 1.5) +
#   geom_point(color = "#e6194b", size = 4) +
#   labs(title = "Cluster Reduction", y = "Number of Clusters", x = "Round") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold", size = 14))
# 
# gridExtra::grid.arrange(p1, p2, ncol = 2)

## ----sankey_basic, fig.height=7, eval=FALSE-----------------------------------
# if (requireNamespace("ggalluvial", quietly = TRUE)) {
#   plot_cluster_sankey(
#     labels_from = labels_over,
#     labels_to = as.character(optim_result$final_labels)
#   )
# }

## ----sankey_custom, fig.height=7, eval=FALSE----------------------------------
# if (requireNamespace("ggalluvial", quietly = TRUE)) {
#   custom_colors <- c(
#     "Type_A1" = "#e6194b", "Type_A2" = "#f58231",
#     "Type_B1" = "#3cb44b", "Type_B2" = "#46f0f0",
#     "Type_C" = "#4363d8", "Type_D" = "#911eb4", "Type_E" = "#f032e6",
#     "1" = "#808080", "2" = "#808080", "3" = "#808080",
#     "4" = "#808080", "5" = "#808080"
#   )
# 
#   plot_cluster_sankey(
#     labels_from = labels_over,
#     labels_to = as.character(optim_result$final_labels),
#     title = "Cluster Merging Flow",
#     colors = custom_colors,
#     alpha = 0.7
#   )
# }

## ----publication_fig, fig.width=14, fig.height=10-----------------------------
# # Create comprehensive figure
# library(gridExtra)
# 
# # Panel A: ROC curves
# p_roc <- plot_roc(result, plot_type = "roc", show_auc = FALSE, legend_position = "none") +
#   labs(title = "A. ROC Curves") +
#   theme(plot.title = element_text(face = "bold", size = 12))
# 
# # Panel B: Confusion heatmap
# p_conf <- plot_confusion_heatmap(result, normalized = "R1", show_values = TRUE, text_size = 3) +
#   labs(title = "B. R1-Normalized Confusion") +
#   theme(plot.title = element_text(face = "bold", size = 12))
# 
# # Panel C: Per-cluster accuracy
# acc_df <- data.frame(
#   Cluster = factor(names(result$per_class_accuracy),
#                    levels = names(sort(result$per_class_accuracy))),
#   Accuracy = result$per_class_accuracy
# )
# p_acc <- ggplot(acc_df, aes(x = Cluster, y = Accuracy)) +
#   geom_col(fill = "#3cb44b", width = 0.6) +
#   geom_hline(yintercept = result$accuracy, linetype = "dashed", color = "red") +
#   coord_flip() +
#   labs(title = "C. Per-Cluster Accuracy", y = "Accuracy", x = NULL) +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold", size = 12)) +
#   ylim(0, 1)
# 
# # Panel D: Metrics summary
# metrics <- data.frame(
#   Metric = c("Test Accuracy", "CV Accuracy", "Max R1", "Max R2"),
#   Value = c(result$accuracy, result$cv_accuracy, result$max_r1, result$max_r2)
# )
# p_metrics <- ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
#   geom_col(width = 0.6) +
#   geom_text(aes(label = sprintf("%.3f", Value)), vjust = -0.3) +
#   scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#d62728", "#2ca02c")) +
#   labs(title = "D. Assessment Metrics", y = "Value", x = NULL) +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold", size = 12),
#         legend.position = "none",
#         axis.text.x = element_text(angle = 30, hjust = 1)) +
#   ylim(0, 1.1)
# 
# # Arrange panels
# grid.arrange(p_roc, p_conf, p_acc, p_metrics,
#              ncol = 2, nrow = 2,
#              top = grid::textGrob("Clustering Assessment Overview",
#                                   gp = grid::gpar(fontface = "bold", fontsize = 16)))

## ----save_example, eval=FALSE-------------------------------------------------
# # Save individual plots
# ggsave("roc_curves.pdf", plot_roc(result), width = 8, height = 6)
# ggsave("confusion_matrix.png", plot_confusion_heatmap(result),
#        width = 7, height = 6, dpi = 300)
# 
# # Save combined figure
# combined_fig <- grid.arrange(p_roc, p_conf, p_acc, p_metrics, ncol = 2)
# ggsave("assessment_overview.pdf", combined_fig, width = 14, height = 10)

## ----custom_theme-------------------------------------------------------------
# # Create a custom theme
# theme_scClustEval <- function() {
#   theme_minimal() +
#     theme(
#       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
#       plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
#       axis.title = element_text(face = "bold", size = 11),
#       axis.text = element_text(size = 10),
#       legend.title = element_text(face = "bold"),
#       panel.grid.minor = element_blank(),
#       strip.text = element_text(face = "bold", size = 11)
#     )
# }
# 
# # Apply custom theme
# plot_roc(result, plot_type = "roc") +
#   theme_scClustEval() +
#   labs(title = "ROC Analysis with Custom Theme")

