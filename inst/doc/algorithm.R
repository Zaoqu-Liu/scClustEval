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

## ----algorithm_flow, echo=FALSE, fig.height=4---------------------------------
# library(ggplot2)
# 
# steps <- data.frame(
#   step = factor(1:5, labels = c(
#     "1. Data Split",
#     "2. Train Classifier",
#     "3. Cross-Validate",
#     "4. Compute Confusion",
#     "5. Merge Clusters"
#   )),
#   description = c(
#     "Stratified train/test split",
#     "Multi-class classification",
#     "K-fold CV on training set",
#     "Build confusion matrix",
#     "Graph-based merging"
#   ),
#   x = 1:5,
#   y = rep(1, 5)
# )
# 
# ggplot(steps, aes(x = x, y = y)) +
#   geom_segment(aes(xend = x, yend = 0), color = "#3cb44b", size = 2) +
#   geom_point(size = 8, color = "#3cb44b") +
#   geom_text(aes(label = step), y = -0.3, size = 3.5, fontface = "bold") +
#   geom_text(aes(label = description), y = 0.3, size = 3) +
#   theme_void() +
#   ylim(-0.8, 0.8) +
#   labs(title = "Self-Projection Algorithm Pipeline") +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

## ----r1_example, echo=TRUE----------------------------------------------------
# library(scClustEval)
# 
# # Example confusion matrix
# cmat <- matrix(c(
#   90, 5, 3,    # Cluster A: 90 correct, 5 misclassified as B, 3 as C
#   8, 85, 7,   # Cluster B: 8 misclassified as A, 85 correct, 7 as C
#   2, 12, 86   # Cluster C: 2 as A, 12 as B, 86 correct
# ), nrow = 3, byrow = TRUE)
# rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")
# 
# # Compute R1 normalization
# r1_mat <- normalize_confmat_r1(cmat)
# print(round(r1_mat, 3))

## ----r2_example, echo=TRUE----------------------------------------------------
# # Compute R2 normalization
# r2_mat <- normalize_confmat_r2(cmat)
# print(round(r2_mat, 4))

## ----confusion_viz, echo=TRUE, fig.height=5, fig.width=10---------------------
# library(ggplot2)
# library(gridExtra)
# 
# # Create data for visualization
# df_r1 <- as.data.frame(as.table(r1_mat))
# colnames(df_r1) <- c("True", "Predicted", "R1")
# 
# df_r2 <- as.data.frame(as.table(r2_mat))
# colnames(df_r2) <- c("True", "Predicted", "R2")
# 
# p1 <- ggplot(df_r1, aes(x = Predicted, y = True, fill = R1)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = sprintf("%.2f", R1)), color = "white", size = 5) +
#   scale_fill_gradient(low = "gray90", high = "#d62728") +
#   labs(title = "R1 Normalization", subtitle = "Pairwise confusion rate") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))
# 
# p2 <- ggplot(df_r2, aes(x = Predicted, y = True, fill = R2)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = sprintf("%.3f", R2)), color = "white", size = 5) +
#   scale_fill_gradient(low = "gray90", high = "#1f77b4") +
#   labs(title = "R2 Normalization", subtitle = "Global confusion proportion") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))
# 
# grid.arrange(p1, p2, ncol = 2)

## ----merging_demo, echo=TRUE, fig.height=4------------------------------------
# # Demonstrate adjacency matrix clustering
# adj <- matrix(c(
#   0, 0.3, 0.02, 0.01,
#   0.3, 0, 0.25, 0.02,
#   0.02, 0.25, 0, 0.01,
#   0.01, 0.02, 0.01, 0
# ), nrow = 4, byrow = TRUE)
# rownames(adj) <- colnames(adj) <- c("C1", "C2", "C3", "C4")
# 
# # Cluster at threshold 0.1
# groups <- cluster_adjacency_matrix(adj, cutoff = 0.1, resolution = 1.0)
# names(groups) <- c("C1", "C2", "C3", "C4")
# print(groups)
# 
# # Visualize
# adj_df <- as.data.frame(as.table(adj))
# colnames(adj_df) <- c("From", "To", "Weight")
# 
# ggplot(adj_df, aes(x = To, y = From, fill = Weight)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = sprintf("%.2f", Weight)), size = 4) +
#   scale_fill_gradient2(low = "white", mid = "#ffbb78", high = "#d62728", midpoint = 0.15) +
#   geom_hline(yintercept = 2.5, linetype = "dashed", color = "black", size = 1) +
#   geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", size = 1) +
#   labs(title = "Confusion-Based Adjacency Matrix",
#        subtitle = "Dashed lines show detected communities (C1,C2 | C3,C4)",
#        fill = "R1 Value") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))

## ----optimization_flow, echo=FALSE, fig.height=5------------------------------
# # Simulated optimization trajectory
# rounds <- 1:8
# accuracy <- c(0.72, 0.78, 0.83, 0.86, 0.88, 0.91, 0.93, 0.94)
# clusters <- c(15, 12, 10, 8, 7, 6, 5, 5)
# 
# df <- data.frame(
#   Round = rounds,
#   Accuracy = accuracy,
#   Clusters = clusters
# )
# 
# p1 <- ggplot(df, aes(x = Round, y = Accuracy)) +
#   geom_line(color = "#3cb44b", size = 1.5) +
#   geom_point(color = "#3cb44b", size = 3) +
#   geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
#   annotate("text", x = 7, y = 0.92, label = "Target: 90%", color = "red") +
#   labs(title = "Accuracy Progression", y = "Test Accuracy") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))
# 
# p2 <- ggplot(df, aes(x = Round, y = Clusters)) +
#   geom_line(color = "#e6194b", size = 1.5) +
#   geom_point(color = "#e6194b", size = 3) +
#   labs(title = "Cluster Count", y = "Number of Clusters") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))
# 
# grid.arrange(p1, p2, ncol = 2)

