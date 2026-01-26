// confusion_matrix.cpp
// Fast confusion matrix computation and normalization
// Author: Zaoqu Liu

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Compute confusion matrix (C++ version)
//'
//' @param y_true True labels (0-indexed integers)
//' @param y_pred Predicted labels (0-indexed integers)
//' @param n_classes Number of classes
//' @return Confusion matrix (rows = true, cols = predicted)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::imat confusion_matrix_cpp(const arma::ivec& y_true,
                                const arma::ivec& y_pred,
                                int n_classes) {
    arma::imat cmat(n_classes, n_classes, fill::zeros);

    int n = y_true.n_elem;
    for (int i = 0; i < n; i++) {
        int true_label = y_true(i);
        int pred_label = y_pred(i);

        if (true_label >= 0 && true_label < n_classes &&
            pred_label >= 0 && pred_label < n_classes) {
            cmat(true_label, pred_label)++;
        }
    }

    return cmat;
}

//' R1 normalize confusion matrix (C++ version)
//'
//' Computes the maximum confusion ratio for each cluster pair.
//'
//' @param cmat Confusion matrix
//' @return R1-normalized confusion matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat normalize_confmat_r1_cpp(const arma::mat& cmat) {
    int n = cmat.n_rows;
    arma::mat xmat(n, n, fill::zeros);

    // Get diagonal values
    arma::vec diag_vals = cmat.diag();

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double val1 = 0.0, val2 = 0.0;

            if (diag_vals(j) > 0) {
                val1 = cmat(i, j) / diag_vals(j);
            }
            if (diag_vals(i) > 0) {
                val2 = cmat(j, i) / diag_vals(i);
            }

            double max_val = std::max(val1, val2);
            xmat(i, j) = max_val;
            xmat(j, i) = max_val;
        }
    }

    return xmat;
}

//' R2 normalize confusion matrix (C++ version)
//'
//' Computes the sum of bidirectional confusions divided by total count.
//'
//' @param cmat Confusion matrix
//' @return R2-normalized confusion matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat normalize_confmat_r2_cpp(const arma::mat& cmat) {
    int n = cmat.n_rows;
    double total = accu(cmat);

    if (total <= 0) {
        return arma::mat(n, n, fill::zeros);
    }

    arma::mat xmat = (cmat + cmat.t()) / total;

    // Set diagonal to zero
    xmat.diag().zeros();

    return xmat;
}

//' Compute per-cluster accuracy (C++ version)
//'
//' @param cmat Confusion matrix
//' @return Vector of per-cluster accuracies
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec per_cluster_accuracy_cpp(const arma::mat& cmat) {
    int n = cmat.n_rows;
    arma::vec accuracy(n);

    arma::vec diag_vals = cmat.diag();
    arma::vec row_sums = sum(cmat, 1);

    for (int i = 0; i < n; i++) {
        if (row_sums(i) > 0) {
            accuracy(i) = diag_vals(i) / row_sums(i);
        } else {
            accuracy(i) = 0.0;
        }
    }

    return accuracy;
}

//' Aggregate confusion matrices by taking element-wise minimum
//'
//' @param mat_list List of confusion matrices
//' @return Aggregated matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat aggregate_confmat_min_cpp(const Rcpp::List& mat_list) {
    if (mat_list.size() == 0) {
        Rcpp::stop("Empty matrix list");
    }

    arma::mat result = Rcpp::as<arma::mat>(mat_list[0]);

    for (int k = 1; k < mat_list.size(); k++) {
        arma::mat current = Rcpp::as<arma::mat>(mat_list[k]);
        result = arma::min(result, current);
    }

    return result;
}

//' Find pairs exceeding threshold in symmetric matrix
//'
//' @param mat Symmetric matrix
//' @param threshold Threshold value
//' @return Matrix with 2 columns (i, j) for pairs exceeding threshold
//'
//' @keywords internal
// [[Rcpp::export]]
arma::imat find_pairs_above_threshold_cpp(const arma::mat& mat, double threshold) {
    int n = mat.n_rows;
    std::vector<int> row_idx, col_idx;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (mat(i, j) > threshold) {
                row_idx.push_back(i);
                col_idx.push_back(j);
            }
        }
    }

    int n_pairs = row_idx.size();
    arma::imat result(n_pairs, 2);

    for (int k = 0; k < n_pairs; k++) {
        result(k, 0) = row_idx[k];
        result(k, 1) = col_idx[k];
    }

    return result;
}
