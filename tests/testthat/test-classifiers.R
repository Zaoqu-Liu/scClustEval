# Tests for classifier functions

test_that("create_classifier creates valid LR classifier", {
    skip_if_not_installed("glmnet")

    clf <- create_classifier("LR", penalty = "l1")

    expect_s3_class(clf, "scClustEval_classifier")
    expect_equal(clf$type, "LR")

    # Test fit and predict
    set.seed(42)
    X <- matrix(rnorm(100 * 10), nrow = 100)
    y <- factor(rep(c("A", "B"), each = 50))

    clf$fit(X, y)

    pred <- clf$predict(X)
    expect_equal(length(pred), 100)
    expect_true(all(pred %in% c("A", "B")))

    probs <- clf$predict_prob(X)
    expect_equal(dim(probs), c(100, 2))
    expect_true(all(probs >= 0 & probs <= 1))
    # Probabilities should sum to 1 for each row
    expect_true(all(abs(rowSums(probs) - 1) < 1e-6))
})

test_that("create_classifier creates valid RF classifier", {
    skip_if_not_installed("randomForest")

    clf <- create_classifier("RF", n_trees = 50)

    expect_s3_class(clf, "scClustEval_classifier")
    expect_equal(clf$type, "RF")

    set.seed(42)
    X <- matrix(rnorm(100 * 10), nrow = 100)
    y <- factor(rep(c("A", "B"), each = 50))

    clf$fit(X, y)

    pred <- clf$predict(X)
    expect_equal(length(pred), 100)

    probs <- clf$predict_prob(X)
    expect_equal(dim(probs), c(100, 2))
})

test_that("create_classifier creates valid DT classifier", {
    skip_if_not_installed("rpart")

    clf <- create_classifier("DT")

    expect_s3_class(clf, "scClustEval_classifier")
    expect_equal(clf$type, "DT")

    set.seed(42)
    X <- matrix(rnorm(100 * 10), nrow = 100)
    y <- factor(rep(c("A", "B"), each = 50))

    clf$fit(X, y)

    pred <- clf$predict(X)
    expect_equal(length(pred), 100)
})

test_that("get_available_classifiers returns correct info", {
    info <- get_available_classifiers()

    expect_s3_class(info, "data.frame")
    expect_true("name" %in% colnames(info))
    expect_true("package" %in% colnames(info))
    expect_true("installed" %in% colnames(info))

    expect_true("LR" %in% info$name)
    expect_true("RF" %in% info$name)
    expect_true("SVM" %in% info$name)
})

test_that("classifier handles multiclass correctly", {
    skip_if_not_installed("glmnet")

    set.seed(42)
    n_per_class <- 30
    X1 <- matrix(rnorm(n_per_class * 10, mean = 0), nrow = n_per_class)
    X2 <- matrix(rnorm(n_per_class * 10, mean = 3), nrow = n_per_class)
    X3 <- matrix(rnorm(n_per_class * 10, mean = 6), nrow = n_per_class)
    X <- rbind(X1, X2, X3)
    y <- factor(c(rep("A", n_per_class), rep("B", n_per_class), rep("C", n_per_class)))

    clf <- create_classifier("LR", penalty = "l2")
    clf$fit(X, y)

    pred <- clf$predict(X)
    expect_true(all(pred %in% c("A", "B", "C")))

    probs <- clf$predict_prob(X)
    expect_equal(ncol(probs), 3)
    expect_true(all(colnames(probs) %in% c("A", "B", "C")))
})

test_that("classifier handles sparse matrices", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("Matrix")

    set.seed(42)
    X_dense <- matrix(rnorm(100 * 10), nrow = 100)
    X_sparse <- Matrix::Matrix(X_dense, sparse = TRUE)
    y <- factor(rep(c("A", "B"), each = 50))

    clf <- create_classifier("LR")
    clf$fit(X_sparse, y)

    pred <- clf$predict(X_sparse)
    expect_equal(length(pred), 100)

    probs <- clf$predict_prob(X_sparse)
    expect_equal(dim(probs), c(100, 2))
})
