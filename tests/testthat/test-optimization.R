# Tests for optimization functions

test_that("sc_optimize works with simple data", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("igraph")

    set.seed(42)

    # Create data where clusters A and B are similar, C is different
    n_per_class <- 30
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 0.5), nrow = n_per_class)  # Similar to X1
    X3 <- matrix(rnorm(n_per_class * 10, mean = 10), nrow = n_per_class)   # Very different

    X <- rbind(X1, X2, X3)
    labels <- factor(c(rep("A", n_per_class), rep("B", n_per_class), rep("C", n_per_class)))

    result <- sc_optimize(
        X = X,
        labels = labels,
        classifier = "LR",
        n_iter = 1,
        cv = 0,
        r1_cutoff = 0.3,
        verbose = FALSE
    )

    expect_true("new_labels" %in% names(result))
    expect_true("old_labels" %in% names(result))
    expect_true("accuracy" %in% names(result))
    expect_true("converged" %in% names(result))

    # Number of clusters should decrease or stay the same
    expect_true(result$n_clusters_after <= result$n_clusters_before)
})

test_that("sc_optimize_all converges", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("igraph")

    set.seed(123)

    # Create well-separated data
    n_per_class <- 40
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 5), nrow = n_per_class)
    X3 <- matrix(rnorm(n_per_class * 10, mean = 10), nrow = n_per_class)

    X <- rbind(X1, X2, X3)
    labels <- factor(c(rep("A", n_per_class), rep("B", n_per_class), rep("C", n_per_class)))

    result <- sc_optimize_all(
        X = X,
        labels = labels,
        min_accuracy = 0.8,
        max_rounds = 3,
        classifier = "LR",
        cv = 0,
        n_iter = 1,
        verbose = FALSE
    )

    expect_s3_class(result, "scClustEval_optim")
    expect_true("final_labels" %in% names(result))
    expect_true("accuracy_history" %in% names(result))
    expect_true("final_accuracy" %in% names(result))

    # Final accuracy should be >= min_accuracy or convergence reached
    expect_true(result$final_accuracy >= 0.7 || result$converged)
})

test_that("sc_optimize preserves cell count", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("igraph")

    set.seed(42)
    n_cells <- 100
    X <- matrix(rnorm(n_cells * 10), nrow = n_cells)
    labels <- factor(sample(c("A", "B", "C", "D"), n_cells, replace = TRUE))

    result <- sc_optimize(
        X = X,
        labels = labels,
        classifier = "LR",
        n_iter = 1,
        cv = 0,
        verbose = FALSE
    )

    # Cell count should be preserved
    expect_equal(length(result$new_labels), n_cells)
    expect_equal(length(result$old_labels), n_cells)
})

test_that("optimization history tracking works", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("igraph")

    set.seed(42)
    n_per_class <- 30
    X <- matrix(rnorm(n_per_class * 3 * 10), nrow = n_per_class * 3)
    labels <- factor(rep(c("A", "B", "C"), each = n_per_class))

    result <- sc_optimize_all(
        X = X,
        labels = labels,
        min_accuracy = 0.99,  # Hard to reach
        max_rounds = 2,
        classifier = "LR",
        cv = 0,
        n_iter = 1,
        verbose = FALSE
    )

    # Should have history entries
    expect_true(length(result$accuracy_history) >= 1)
    expect_true(length(result$n_clusters_history) >= 2)  # Includes initial

    # History should be consistent
    expect_equal(result$total_rounds, length(result$accuracy_history))
})
