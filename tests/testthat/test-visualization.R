# Tests for visualization functions

test_that("plot_roc returns ggplot object", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("pROC")

    set.seed(42)
    n_per_class <- 50
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 5), nrow = n_per_class)
    X <- rbind(X1, X2)
    y <- factor(c(rep("A", n_per_class), rep("B", n_per_class)))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    p <- plot_roc(result)
    expect_s3_class(p, "gg")
})

test_that("plot_confusion_heatmap works with result object", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 50
    X <- rbind(
        matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class),
        matrix(rnorm(n_per_class * 10, mean = 5), nrow = n_per_class)
    )
    y <- factor(c(rep("A", n_per_class), rep("B", n_per_class)))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    # Test all normalization modes
    p1 <- plot_confusion_heatmap(result, normalized = "raw")
    p2 <- plot_confusion_heatmap(result, normalized = "R1")
    p3 <- plot_confusion_heatmap(result, normalized = "R2")

    expect_s3_class(p1, "gg")
    expect_s3_class(p2, "gg")
    expect_s3_class(p3, "gg")
})

test_that("plot_confusion_heatmap works with matrix directly", {
    skip_if_not_installed("ggplot2")

    cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")

    p <- plot_confusion_heatmap(cmat)
    expect_s3_class(p, "gg")
})

test_that("plot_optimization_history works", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 40
    X <- rbind(
        matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class),
        matrix(rnorm(n_per_class * 10, mean = 5), nrow = n_per_class),
        matrix(rnorm(n_per_class * 10, mean = 10), nrow = n_per_class)
    )
    y <- factor(rep(c("A", "B", "C"), each = n_per_class))

    result <- sc_optimize_all(
        X = X,
        labels = y,
        min_accuracy = 0.99,
        max_rounds = 2,
        classifier = "LR",
        cv = 0,
        n_iter = 1,
        verbose = FALSE
    )

    p <- plot_optimization_history(result)
    expect_s3_class(p, "gg")
})

test_that("plot.scClustEval method works", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n <- 100
    X <- matrix(rnorm(n * 10), nrow = n)
    y <- factor(rep(c("A", "B"), each = 50))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    p1 <- plot(result, type = "roc")
    p2 <- plot(result, type = "confusion")

    expect_s3_class(p1, "gg")
    expect_s3_class(p2, "gg")
})

test_that("print and summary methods work", {
    skip_if_not_installed("glmnet")

    set.seed(42)
    n <- 100
    X <- matrix(rnorm(n * 10), nrow = n)
    y <- factor(rep(c("A", "B"), each = 50))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    # Should not error
    expect_output(print(result))
    expect_output(summary(result))
})
