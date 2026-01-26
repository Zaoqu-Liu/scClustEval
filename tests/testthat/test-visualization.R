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

test_that("plot_cluster_links returns ggplot object", {
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 30
    n_clusters <- 3

    # Create embeddings with distinct cluster centers
    embeddings <- matrix(rnorm(n_per_class * n_clusters * 2), ncol = 2)
    labels <- factor(rep(1:n_clusters, each = n_per_class))

    # Separate cluster centers
    for (i in 1:n_clusters) {
        idx <- which(labels == i)
        embeddings[idx, 1] <- embeddings[idx, 1] + (i - 1) * 5
    }

    # Create a fake confusion matrix
    cmat <- matrix(c(0, 0.3, 0.1, 0.3, 0, 0.2, 0.1, 0.2, 0), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- levels(labels)

    p <- plot_cluster_links(embeddings, labels, cmat, threshold = 0.1)
    expect_s3_class(p, "gg")
})

test_that("plot_cluster_centers returns ggplot object", {
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 30
    n_clusters <- 3

    embeddings <- matrix(rnorm(n_per_class * n_clusters * 2), ncol = 2)
    labels <- factor(rep(1:n_clusters, each = n_per_class))

    p <- plot_cluster_centers(embeddings, labels)
    expect_s3_class(p, "gg")
})

test_that("plot_embedding_with_links returns ggplot object", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 30
    n_clusters <- 3

    # Create separable data
    X <- rbind(
        matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class),
        matrix(rnorm(n_per_class * 10, mean = 3), nrow = n_per_class),
        matrix(rnorm(n_per_class * 10, mean = 6), nrow = n_per_class)
    )
    labels <- factor(rep(1:n_clusters, each = n_per_class))

    # Create embeddings
    embeddings <- matrix(rnorm(n_per_class * n_clusters * 2), ncol = 2)
    for (i in 1:n_clusters) {
        idx <- which(labels == i)
        embeddings[idx, 1] <- embeddings[idx, 1] + (i - 1) * 5
    }

    # Run assessment
    result <- sc_assessment(X, labels, classifier = "LR", cv = 0, verbose = FALSE)

    # Test combined plot
    p <- plot_embedding_with_links(
        embeddings, labels, result$r1_normalized,
        threshold = 0.01
    )
    expect_s3_class(p, "gg")
})

test_that("plot_cluster_links works with add_to_plot parameter", {
    skip_if_not_installed("ggplot2")

    set.seed(42)
    n_per_class <- 30
    n_clusters <- 3

    embeddings <- matrix(rnorm(n_per_class * n_clusters * 2), ncol = 2)
    labels <- factor(rep(1:n_clusters, each = n_per_class))

    cmat <- matrix(c(0, 0.3, 0.1, 0.3, 0, 0.2, 0.1, 0.2, 0), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- levels(labels)

    # Create base plot
    base_plot <- ggplot2::ggplot() + ggplot2::theme_minimal()

    # Add links to existing plot
    p <- plot_cluster_links(embeddings, labels, cmat,
                            threshold = 0.05, add_to_plot = base_plot)
    expect_s3_class(p, "gg")
})
