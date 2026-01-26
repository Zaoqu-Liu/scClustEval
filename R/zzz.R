#' @title Package Startup Functions
#' @description Functions called when the package is loaded
#' @name zzz
NULL

# Package-level environment for internal state
.scClustEval_env <- new.env(parent = emptyenv())

# Default color palette (20 colors, based on Tableau palette)
.color_palette_20 <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
)

# Extended color palette (100+ colors for large datasets)
.color_palette_extended <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
    "#5254a3", "#8ca252", "#bd9e39", "#ad494a", "#a55194",
    "#6b6ecf", "#b5cf6b", "#e7ba52", "#d6616b", "#ce6dbd",
    "#9c9ede", "#cedb9c", "#e7cb94", "#e7969c", "#de9ed6",
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#e6beff", "#aa6e28", "#800000", "#aaffc3",
    "#808000", "#ffd8b1", "#000080", "#808080", "#000000",
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
    "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6",
    "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601",
    "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900",
    "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
    "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
    "#300018", "#0AA6D8", "#013349", "#00846F", "#372101"
)

.onLoad <- function(libname, pkgname) {
    # Set default options
    op <- options()
    op.scClustEval <- list(
        scClustEval.verbose = TRUE,
        scClustEval.n_cores = 1L
    )

    toset <- !(names(op.scClustEval) %in% names(op))
    if (any(toset)) options(op.scClustEval[toset])

    # Initialize package environment
    .scClustEval_env$version <- utils::packageVersion(pkgname)

    invisible()
}

.onAttach <- function(libname, pkgname) {
    version <- utils::packageVersion(pkgname)
    packageStartupMessage(sprintf(
        "scClustEval v%s - Single Cell Clustering Evaluation Framework",
        version
    ))
    packageStartupMessage("Author: Zaoqu Liu <liuzaoqu@163.com>")
    packageStartupMessage("GitHub: https://github.com/Zaoqu-Liu/scClustEval")
}

.onUnload <- function(libpath) {
    # Clean up
    tryCatch({
        library.dynam.unload("scClustEval", libpath)
    }, error = function(e) NULL)
}
