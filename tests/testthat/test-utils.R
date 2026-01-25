# Tests for utility functions

test_that("make_unique_names works correctly", {
    # Basic case
    result <- make_unique_names(c("A", "B", "A", "C", "A"))
    expect_equal(result, c("A", "B", "A_1", "C", "A_2"))

    # No duplicates
    result <- make_unique_names(c("A", "B", "C"))
    expect_equal(result, c("A", "B", "C"))

    # All same
    result <- make_unique_names(c("X", "X", "X"))
    expect_equal(result, c("X", "X_1", "X_2"))
})

test_that("train_test_split_stratified works correctly", {
    set.seed(42)
    X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
    y <- factor(rep(c("A", "B"), each = 50))

    result <- train_test_split_stratified(X, y, test_size = 0.5, seed = 42)

    # Check structure
    expect_true(all(c("X_train", "X_test", "y_train", "y_test") %in% names(result)))

    # Check dimensions
    expect_equal(nrow(result$X_train) + nrow(result$X_test), nrow(X))
    expect_equal(length(result$y_train) + length(result$y_test), length(y))

    # Check stratification (approximately equal split per class)
    train_table <- table(result$y_train)
    expect_true(abs(train_table["A"] - train_table["B"]) <= 5)
})
