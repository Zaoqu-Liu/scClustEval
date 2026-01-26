// distance.cpp
// Fast distance matrix computation using RcppArmadillo
// Author: Zaoqu Liu
// Note: No OpenMP to avoid conflicts with other packages

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Euclidean distance matrix (C++ version)
//'
//' @param X Numeric matrix (rows = samples, cols = features)
//' @return Distance matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat euclidean_distance_cpp(const arma::mat& X) {
    int n = X.n_rows;
    arma::mat D(n, n, fill::zeros);

    // Compute squared norms
    arma::vec sq_norms = sum(X % X, 1);

    // D_ij = ||x_i||^2 + ||x_j||^2 - 2 * x_i . x_j
    arma::mat XXT = X * X.t();

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d_sq = sq_norms(i) + sq_norms(j) - 2.0 * XXT(i, j);
            // Handle numerical issues
            d_sq = std::max(0.0, d_sq);
            D(i, j) = std::sqrt(d_sq);
            D(j, i) = D(i, j);
        }
    }

    return D;
}

//' Compute pairwise distances between two matrices (C++ version)
//'
//' @param X First matrix (n x p)
//' @param Y Second matrix (m x p)
//' @return Distance matrix (n x m)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat pairwise_distance_cpp(const arma::mat& X, const arma::mat& Y) {
    int n = X.n_rows;
    int m = Y.n_rows;
    arma::mat D(n, m);

    arma::vec X_sq = sum(X % X, 1);
    arma::vec Y_sq = sum(Y % Y, 1);
    arma::mat XY = X * Y.t();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double d_sq = X_sq(i) + Y_sq(j) - 2.0 * XY(i, j);
            d_sq = std::max(0.0, d_sq);
            D(i, j) = std::sqrt(d_sq);
        }
    }

    return D;
}

//' Compute cluster centroids (C++ version)
//'
//' @param X Expression matrix (cells x features)
//' @param labels Cluster labels (0-indexed integers)
//' @param n_clusters Number of clusters
//' @return Centroid matrix (clusters x features)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_centroids_cpp(const arma::mat& X,
                                const arma::ivec& labels,
                                int n_clusters) {
    int n_features = X.n_cols;
    arma::mat centroids(n_clusters, n_features, fill::zeros);
    arma::ivec counts(n_clusters, fill::zeros);

    int n_cells = X.n_rows;

    // Sum up features for each cluster
    for (int i = 0; i < n_cells; i++) {
        int cl = labels(i);
        if (cl >= 0 && cl < n_clusters) {
            centroids.row(cl) += X.row(i);
            counts(cl)++;
        }
    }

    // Compute means
    for (int c = 0; c < n_clusters; c++) {
        if (counts(c) > 0) {
            centroids.row(c) /= counts(c);
        }
    }

    return centroids;
}

//' Compute cosine similarity matrix (C++ version)
//'
//' @param X Numeric matrix (rows = samples, cols = features)
//' @return Similarity matrix with values between -1 and 1
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat cosine_similarity_cpp(const arma::mat& X) {
    // Normalize rows
    arma::mat X_norm = X;
    arma::vec row_norms = sqrt(sum(X % X, 1));

    for (arma::uword i = 0; i < X.n_rows; i++) {
        if (row_norms(i) > 0) {
            X_norm.row(i) /= row_norms(i);
        }
    }

    // Cosine similarity = normalized dot product
    return X_norm * X_norm.t();
}
