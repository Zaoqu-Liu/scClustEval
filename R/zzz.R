#' @title Package Startup Functions
#' @description Functions called when the package is loaded
#' @name zzz
NULL

.onLoad <- function(libname, pkgname) {
    # Set default options
    op <- options()
    op.scClustEval <- list(
        scClustEval.verbose = TRUE,
        scClustEval.n_cores = 1L,
        scClustEval.use_cpp = TRUE
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
    library.dynam.unload("scClustEval", libpath)
}
