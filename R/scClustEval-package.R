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
#' @importFrom Rcpp sourceCpp
#'
"_PACKAGE"

# Package-level variables
.scClustEval_env <- new.env(parent = emptyenv())

#' @title Available color palettes
#' @description Extended color palette for visualization
#' @keywords internal
.color_palette_extended <- c(
    '#e6194b', '#3cb44b', '#ffe119', '#0082c8', '#f58231', '#911eb4',
    '#46f0f0', '#f032e6', '#d2f53c', '#fabebe', '#008080', '#e6beff',
    '#aa6e28', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000080',
    '#808080', '#000000', "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46",
    "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6",
    "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
    "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53",
    "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
    "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B",
    "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
    "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99",
    "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75",
    "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A",
    "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C", "#34362D",
    "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
    "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700",
    "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9",
    "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79",
    "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465",
    "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98",
    "#A4E804", "#324E72"
)

#' @title Default 20-color palette
#' @description Default palette for up to 20 clusters
#' @keywords internal
.color_palette_20 <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
)
