#' @title Classifier Functions for scClustEval
#' @description Unified interface for multiple machine learning classifiers
#' @name classifiers
NULL

#' @title Get available classifiers
#' @description List all available classifiers and their requirements
#'
#' @return A data.frame with classifier information
#'
#' @examples
#' get_available_classifiers()
#'
#' @export
get_available_classifiers <- function() {
    classifiers <- data.frame(
        name = c("LR", "RF", "SVM", "NB", "DT", "XGB", "RANGER"),
        full_name = c(
            "Logistic Regression (glmnet)",
            "Random Forest (randomForest)",
            "Support Vector Machine (e1071)",
            "Naive Bayes (e1071)",
            "Decision Tree (rpart)",
            "XGBoost (xgboost)",
            "Random Forest (ranger)"
        ),
        package = c("glmnet", "randomForest", "e1071", "e1071", "rpart", "xgboost", "ranger"),
        installed = c(
            requireNamespace("glmnet", quietly = TRUE),
            requireNamespace("randomForest", quietly = TRUE),
            requireNamespace("e1071", quietly = TRUE),
            requireNamespace("e1071", quietly = TRUE),
            requireNamespace("rpart", quietly = TRUE),
            requireNamespace("xgboost", quietly = TRUE),
            requireNamespace("ranger", quietly = TRUE)
        ),
        supports_probability = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    return(classifiers)
}

#' @title Create a classifier
#' @description Create a unified classifier object that wraps various ML algorithms
#'
#' @param type Classifier type: "LR", "RF", "SVM", "NB", "DT", "XGB", "RANGER"
#' @param penalty For LR: "l1" (lasso), "l2" (ridge), or "elasticnet"
#' @param alpha For LR: elasticnet mixing parameter (1=lasso, 0=ridge). Default: 1 for L1
#' @param lambda For LR: regularization strength (smaller = more regularization).
#'   If NULL, uses cross-validation to select
#' @param n_trees For RF/RANGER/XGB: number of trees
#' @param max_depth For DT/XGB: maximum tree depth
#' @param kernel For SVM: kernel type ("linear", "radial", "polynomial")
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments passed to the underlying classifier
#'
#' @return A classifier object with fit, predict, and predict_prob methods
#'
#' @details
#' The returned classifier object has the following methods:
#' \describe{
#'   \item{$fit(X, y)}{Train the classifier on data}
#'   \item{$predict(X)}{Get class predictions}
#'   \item{$predict_prob(X)}{Get class probabilities}
#'   \item{$get_coefficients()}{Get model coefficients (for LR)}
#'   \item{$get_importance()}{Get feature importance (for tree-based)}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a logistic regression classifier
#' clf <- create_classifier("LR", penalty = "l1")
#'
#' # Train
#' clf$fit(X_train, y_train)
#'
#' # Predict
#' predictions <- clf$predict(X_test)
#' probabilities <- clf$predict_prob(X_test)
#' }
#'
#' @export
create_classifier <- function(type = "LR",
                              penalty = "l1",
                              alpha = NULL,
                              lambda = NULL,
                              n_trees = 500,
                              max_depth = NULL,
                              kernel = "radial",
                              seed = NULL,
                              ...) {

    type <- toupper(type)
    valid_types <- c("LR", "RF", "SVM", "NB", "DT", "XGB", "RANGER")

    if (!type %in% valid_types) {
        stop(sprintf("Unknown classifier type: %s. Valid types: %s",
                     type, paste(valid_types, collapse = ", ")))
    }

    # Check package availability
    pkg_map <- c(
        "LR" = "glmnet",
        "RF" = "randomForest",
        "SVM" = "e1071",
        "NB" = "e1071",
        "DT" = "rpart",
        "XGB" = "xgboost",
        "RANGER" = "ranger"
    )

    required_pkg <- pkg_map[type]
    if (!requireNamespace(required_pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is required for classifier '%s'. Please install it.",
                     required_pkg, type))
    }

    # Set alpha based on penalty type
    if (type == "LR" && is.null(alpha)) {
        alpha <- switch(penalty,
                        "l1" = 1,
                        "l2" = 0,
                        "elasticnet" = 0.5,
                        1)
    }

    # Create classifier based on type
    clf <- switch(type,
                  "LR" = .create_lr_classifier(alpha, lambda, seed, ...),
                  "RF" = .create_rf_classifier(n_trees, seed, ...),
                  "SVM" = .create_svm_classifier(kernel, seed, ...),
                  "NB" = .create_nb_classifier(seed, ...),
                  "DT" = .create_dt_classifier(max_depth, seed, ...),
                  "XGB" = .create_xgb_classifier(n_trees, max_depth, seed, ...),
                  "RANGER" = .create_ranger_classifier(n_trees, seed, ...))

    clf$type <- type
    clf$seed <- seed

    class(clf) <- c("scClustEval_classifier", class(clf))
    return(clf)
}

#' @title Create Logistic Regression Classifier
#' @keywords internal
.create_lr_classifier <- function(alpha, lambda, seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            # Convert sparse matrix if needed
            if (inherits(X, "sparseMatrix")) {
                X <- as(X, "dgCMatrix")
            }

            # Fit model with cross-validation if lambda is NULL
            if (is.null(lambda)) {
                model <<- glmnet::cv.glmnet(
                    x = X,
                    y = y,
                    family = "multinomial",
                    alpha = alpha,
                    type.measure = "class",
                    nfolds = 5,
                    ...
                )
            } else {
                model <<- glmnet::glmnet(
                    x = X,
                    y = y,
                    family = "multinomial",
                    alpha = alpha,
                    lambda = lambda,
                    ...
                )
            }
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "sparseMatrix")) {
                X <- as(X, "dgCMatrix")
            }

            s_val <- if (inherits(model, "cv.glmnet")) "lambda.min" else model$lambda[1]
            pred <- predict(model, newx = X, s = s_val, type = "class")
            return(as.character(pred[, 1]))
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "sparseMatrix")) {
                X <- as(X, "dgCMatrix")
            }

            s_val <- if (inherits(model, "cv.glmnet")) "lambda.min" else model$lambda[1]
            probs <- predict(model, newx = X, s = s_val, type = "response")

            # Convert 3D array to matrix
            prob_mat <- probs[, , 1]
            colnames(prob_mat) <- classes
            return(prob_mat)
        },

        get_coefficients = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            s_val <- if (inherits(model, "cv.glmnet")) "lambda.min" else model$lambda[1]
            coefs <- coef(model, s = s_val)
            return(coefs)
        },

        get_importance = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            s_val <- if (inherits(model, "cv.glmnet")) "lambda.min" else model$lambda[1]
            coef_list <- coef(model, s = s_val)
            # Sum absolute values across classes
            coef_mat <- do.call(cbind, lapply(coef_list, function(x) as.vector(x[-1])))
            importance <- rowSums(abs(coef_mat))
            return(importance)
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create Random Forest Classifier
#' @keywords internal
.create_rf_classifier <- function(n_trees, seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            # Convert to dense matrix if sparse
            if (inherits(X, "Matrix")) {
                X <- as.matrix(X)
            }

            model <<- randomForest::randomForest(
                x = X,
                y = y,
                ntree = n_trees,
                importance = TRUE,
                ...
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            pred <- predict(model, newdata = X, type = "response")
            return(as.character(pred))
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            probs <- predict(model, newdata = X, type = "prob")
            return(probs)
        },

        get_coefficients = function() {
            warning("Random Forest does not have coefficients. Use get_importance() instead.")
            return(NULL)
        },

        get_importance = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            return(randomForest::importance(model))
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create SVM Classifier
#' @keywords internal
.create_svm_classifier <- function(kernel, seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            if (inherits(X, "Matrix")) X <- as.matrix(X)

            model <<- e1071::svm(
                x = X,
                y = y,
                kernel = kernel,
                probability = TRUE,
                ...
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            pred <- predict(model, newdata = X)
            return(as.character(pred))
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            pred <- predict(model, newdata = X, probability = TRUE)
            probs <- attr(pred, "probabilities")
            # Reorder columns to match classes order
            probs <- probs[, classes, drop = FALSE]
            return(probs)
        },

        get_coefficients = function() {
            warning("SVM does not have coefficients in the traditional sense.")
            return(NULL)
        },

        get_importance = function() {
            warning("Feature importance not directly available for SVM.")
            return(NULL)
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create Naive Bayes Classifier
#' @keywords internal
.create_nb_classifier <- function(seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            if (inherits(X, "Matrix")) X <- as.matrix(X)

            model <<- e1071::naiveBayes(
                x = X,
                y = y,
                ...
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            pred <- predict(model, newdata = X, type = "class")
            return(as.character(pred))
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            probs <- predict(model, newdata = X, type = "raw")
            return(probs)
        },

        get_coefficients = function() {
            warning("Naive Bayes does not have coefficients.")
            return(NULL)
        },

        get_importance = function() {
            warning("Feature importance not available for Naive Bayes.")
            return(NULL)
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create Decision Tree Classifier
#' @keywords internal
.create_dt_classifier <- function(max_depth, seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            if (inherits(X, "Matrix")) X <- as.matrix(X)

            # Create data frame
            df <- as.data.frame(X)
            df$.y <- y

            control_args <- list(...)
            if (!is.null(max_depth)) {
                control_args$maxdepth <- max_depth
            }

            model <<- rpart::rpart(
                .y ~ .,
                data = df,
                method = "class",
                control = do.call(rpart::rpart.control, control_args)
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            df <- as.data.frame(X)
            pred <- predict(model, newdata = df, type = "class")
            return(as.character(pred))
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            df <- as.data.frame(X)
            probs <- predict(model, newdata = df, type = "prob")
            return(probs)
        },

        get_coefficients = function() {
            warning("Decision Tree does not have coefficients. Use get_importance() instead.")
            return(NULL)
        },

        get_importance = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            return(model$variable.importance)
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create XGBoost Classifier
#' @keywords internal
.create_xgb_classifier <- function(n_trees, max_depth, seed, ...) {
    model <- NULL
    classes <- NULL
    class_map <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)
            class_map <<- setNames(seq_along(classes) - 1, classes)

            # Convert labels to numeric
            y_num <- class_map[as.character(y)]

            if (inherits(X, "Matrix")) {
                X <- as(X, "dgCMatrix")
            }

            params <- list(
                objective = "multi:softprob",
                num_class = length(classes),
                eval_metric = "mlogloss",
                ...
            )

            if (!is.null(max_depth)) {
                params$max_depth <- max_depth
            }

            dtrain <- xgboost::xgb.DMatrix(data = X, label = y_num)

            model <<- xgboost::xgb.train(
                params = params,
                data = dtrain,
                nrounds = n_trees,
                verbose = 0
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")

            if (inherits(X, "Matrix")) {
                X <- as(X, "dgCMatrix")
            }

            dtest <- xgboost::xgb.DMatrix(data = X)
            probs <- predict(model, dtest, reshape = TRUE)
            pred_idx <- apply(probs, 1, which.max)
            return(classes[pred_idx])
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")

            if (inherits(X, "Matrix")) {
                X <- as(X, "dgCMatrix")
            }

            dtest <- xgboost::xgb.DMatrix(data = X)
            probs <- predict(model, dtest, reshape = TRUE)
            colnames(probs) <- classes
            return(probs)
        },

        get_coefficients = function() {
            warning("XGBoost does not have coefficients. Use get_importance() instead.")
            return(NULL)
        },

        get_importance = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            return(xgboost::xgb.importance(model = model))
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Create Ranger (Fast Random Forest) Classifier
#' @keywords internal
.create_ranger_classifier <- function(n_trees, seed, ...) {
    model <- NULL
    classes <- NULL

    list(
        fit = function(X, y) {
            if (!is.null(seed)) set.seed(seed)

            y <- as.factor(y)
            classes <<- levels(y)

            if (inherits(X, "Matrix")) X <- as.matrix(X)

            # Create data frame
            df <- as.data.frame(X)
            df$.y <- y

            model <<- ranger::ranger(
                .y ~ .,
                data = df,
                num.trees = n_trees,
                probability = TRUE,
                importance = "impurity",
                ...
            )
            invisible(NULL)
        },

        predict = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            df <- as.data.frame(X)
            pred <- predict(model, data = df)
            probs <- pred$predictions
            pred_idx <- apply(probs, 1, which.max)
            return(classes[pred_idx])
        },

        predict_prob = function(X) {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            if (inherits(X, "Matrix")) X <- as.matrix(X)
            df <- as.data.frame(X)
            pred <- predict(model, data = df)
            probs <- pred$predictions
            # Ensure column names match classes
            if (is.null(colnames(probs))) {
                colnames(probs) <- classes
            }
            return(probs)
        },

        get_coefficients = function() {
            warning("Ranger does not have coefficients. Use get_importance() instead.")
            return(NULL)
        },

        get_importance = function() {
            if (is.null(model)) stop("Model not trained. Call fit() first.")
            return(model$variable.importance)
        },

        get_classes = function() {
            return(classes)
        },

        get_model = function() {
            return(model)
        }
    )
}

#' @title Print method for classifier
#' @param x Classifier object
#' @param ... Additional arguments
#' @export
#' @method print scClustEval_classifier
print.scClustEval_classifier <- function(x, ...) {
    cat("scClustEval Classifier\n")
    cat("  Type:", x$type, "\n")
    if (!is.null(x$get_classes())) {
        cat("  Classes:", paste(x$get_classes(), collapse = ", "), "\n")
        cat("  Trained: Yes\n")
    } else {
        cat("  Trained: No\n")
    }
    invisible(x)
}
