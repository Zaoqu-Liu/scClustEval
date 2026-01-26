#' @title scClustEval: Single Cell Clustering Evaluation and Optimization Framework
#'
#' @description
#' A comprehensive framework for evaluating and optimizing single-cell RNA-seq
#' clustering results using self-projection machine learning approaches.
#'
#' @details
#' The scClustEval package provides tools for:
#'
#' \itemize{
#'   \item \strong{Clustering Assessment}: Evaluate the quality of cell clustering
#'     using self-projection with various machine learning classifiers
#'   \item \strong{Clustering Optimization}: Iteratively merge poorly discriminated

#'     clusters to achieve robust cell type identification
#'   \item \strong{Visualization}: ROC curves, confusion matrices, Sankey diagrams,
#'     and comprehensive assessment plots
#'   \item \strong{Seurat Integration}: Seamless workflow with Seurat objects
#' }
#'
#' The core algorithm works by:
#' \enumerate{
#'   \item Training a classifier to distinguish between clusters
#'   \item Evaluating prediction accuracy via cross-validation and hold-out testing
#'   \item Identifying cluster pairs that are difficult to discriminate
#'   \item Merging confused clusters and iterating until target accuracy is reached
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{sc_assessment}}}{Core function for clustering assessment}
#'   \item{\code{\link{sc_optimize}}}{Single round of clustering optimization}
#'   \item{\code{\link{sc_optimize_all}}}{Full iterative optimization pipeline}
#'   \item{\code{\link{RunAssessment}}}{Seurat-style assessment function}
#'   \item{\code{\link{RunOptimization}}}{Seurat-style optimization function}
#' }
#'
#' @section Classifiers:
#' The package supports multiple classifiers:
#' \itemize{
#'   \item LR: Logistic Regression (L1/L2 regularization)
#'   \item RF: Random Forest
#'   \item SVM: Support Vector Machine
#'   \item NB: Naive Bayes
#'   \item DT: Decision Tree
#'   \item XGB: XGBoost (if installed)
#' }
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#'
#' @references
#' This package is an R implementation inspired by the SCCAF Python package:
#' \url{https://github.com/SCCAF/sccaf}
#'
#' Miao, Z., et al. (2020). Putative cell type discovery from single-cell gene
#' expression data. Nature Methods.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/Zaoqu-Liu/scClustEval}
#'   \item Report bugs at \url{https://github.com/Zaoqu-Liu/scClustEval/issues}
#' }
#'
#' @docType package
#' @name scClustEval-package
#' @aliases scClustEval
#'
#' @useDynLib scClustEval, .registration = TRUE
#' @import ggplot2
#' @import methods
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @importFrom Matrix Matrix t colMeans rowSums colSums
#' @importFrom stats predict coef model.matrix as.formula median quantile sd var setNames
#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion combn head tail
#' @importFrom methods as
#' @importFrom caret createFolds
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain membership
#' @importFrom pROC roc auc
#' @importFrom rpart rpart rpart.control
#' @importFrom grDevices colorRampPalette
#' @importFrom parallel detectCores
#' @importFrom future plan sequential multisession multicore
#'
"_PACKAGE"

# Note: Package-level variables and color palettes are defined in zzz.R
# to ensure proper initialization order during package loading
