# Tests for assessment functions

test_that("sc_assessment works with simple data", {
    skip_if_not_installed("glmnet")

    set.seed(42)
    # Create well-separated clusters
    n_per_class <- 50
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 5), nrow = n_per_class)
    X <- rbind(X1, X2)
    y <- factor(c(rep("A", n_per_class), rep("B", n_per_class)))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    # Check result structure
    expect_s3_class(result, "scClustEval")
    expect_true("accuracy" %in% names(result))
    expect_true("confusion_matrix" %in% names(result))
    expect_true("r1_normalized" %in% names(result))
    expect_true("r2_normalized" %in% names(result))
    expect_true("per_class_accuracy" %in% names(result))

    # Check accuracy (should be high for well-separated data)
    expect_gt(result$accuracy, 0.7)
})

test_that("sc_assessment handles three or more classes", {
    skip_if_not_installed("glmnet")

    set.seed(123)
    n_per_class <- 30
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 3), nrow = n_per_class)
    X3 <- matrix(rnorm(n_per_class * 10, mean = 6), nrow = n_per_class)
    X <- rbind(X1, X2, X3)
    y <- factor(c(rep("A", n_per_class), rep("B", n_per_class), rep("C", n_per_class)))

    result <- sc_assessment(X, y, classifier = "LR", cv = 0, verbose = FALSE)

    expect_s3_class(result, "scClustEval")
    expect_equal(length(result$classes), 3)
    expect_equal(nrow(result$confusion_matrix), 3)
    expect_equal(ncol(result$confusion_matrix), 3)
})

test_that("confusion_matrix calculation is correct", {
    y_true <- c("A", "A", "B", "B", "C", "C")
    y_pred <- c("A", "B", "B", "B", "C", "A")

    cmat <- calc_confusion_matrix(y_true, y_pred)

    expect_equal(dim(cmat), c(3, 3))
    expect_equal(cmat["A", "A"], 1)
    expect_equal(cmat["A", "B"], 1)
    expect_equal(cmat["B", "B"], 2)
    expect_equal(cmat["C", "C"], 1)
    expect_equal(cmat["C", "A"], 1)
})

test_that("normalize_confmat_r1 produces symmetric matrix", {
    cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")

    r1mat <- normalize_confmat_r1(cmat)

    # Should be symmetric
    expect_equal(r1mat, t(r1mat))

    # Diagonal should be zero (ignore names)
    expect_equal(as.numeric(diag(r1mat)), c(0, 0, 0))
})

test_that("normalize_confmat_r2 produces symmetric matrix with bounded values", {
    cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")

    r2mat <- normalize_confmat_r2(cmat)

    # Should be symmetric
    expect_equal(r2mat, t(r2mat))

    # Diagonal should be zero (ignore names)
    expect_equal(as.numeric(diag(r2mat)), c(0, 0, 0))

    # Values should be bounded [0, 1]
    expect_true(all(r2mat >= 0 & r2mat <= 1))
})

test_that("per_cluster_accuracy computes correctly", {
    cmat <- matrix(c(80, 10, 10, 5, 90, 5, 15, 0, 85), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")

    acc <- per_cluster_accuracy(cmat)

    expect_equal(length(acc), 3)
    # Use unname to ignore the names attribute
    expect_equal(unname(acc["A"]), 80 / 100)
    expect_equal(unname(acc["B"]), 90 / 100)
    expect_equal(unname(acc["C"]), 85 / 100)
})

test_that("train_test_split_stratified maintains class proportions", {
    set.seed(42)
    n <- 100
    X <- matrix(rnorm(n * 5), nrow = n)
    y <- factor(c(rep("A", 40), rep("B", 30), rep("C", 30)))

    split <- train_test_split_stratified(X, y, test_size = 0.3)

    expect_equal(length(split$y_train) + length(split$y_test), n)
    expect_equal(nrow(split$X_train) + nrow(split$X_test), n)

    # Check all classes are present in both splits
    expect_true(all(c("A", "B", "C") %in% split$y_train))
    expect_true(all(c("A", "B", "C") %in% split$y_test))
})

test_that("make_unique_names handles duplicates correctly", {
    x <- c("A", "B", "A", "C", "A", "B")
    result <- make_unique_names(x)

    expect_equal(length(unique(result)), length(x))
    expect_true(all(c("A", "B", "C") %in% result))
})
