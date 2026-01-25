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

    # Check accuracy (should be high for well-separated data)
    expect_gt(result$accuracy, 0.7)
})

test_that("confusion_matrix calculation is correct", {
    y_true <- c("A", "A", "B", "B", "C", "C")
    y_pred <- c("A", "B", "B", "B", "C", "A")

    cmat <- calc_confusion_matrix(y_true, y_pred)

    expect_equal(dim(cmat), c(3, 3))
    expect_equal(cmat["A", "A"], 1)  # 1 correct A prediction
    expect_equal(cmat["A", "B"], 1)  # 1 A misclassified as B
    expect_equal(cmat["B", "B"], 2)  # 2 correct B predictions
    expect_equal(cmat["C", "C"], 1)  # 1 correct C prediction
    expect_equal(cmat["C", "A"], 1)  # 1 C misclassified as A
})

test_that("normalize_confmat_r1 produces symmetric matrix", {
    cmat <- matrix(c(90, 5, 5, 3, 85, 12, 2, 10, 88), nrow = 3)
    rownames(cmat) <- colnames(cmat) <- c("A", "B", "C")

    r1mat <- normalize_confmat_r1(cmat)

    # Should be symmetric
    expect_equal(r1mat, t(r1mat))

    # Diagonal should be zero
    expect_equal(diag(r1mat), c(0, 0, 0))
})
