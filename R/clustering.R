#' @title Clustering Functions for scClustEval
#' @description Functions for graph-based clustering and cluster manipulation
#' @name clustering
NULL

#' @title Cluster adjacency matrix
#' @description Cluster groups based on an adjacency/confusion matrix using Louvain
#'
#' @param adj_matrix Adjacency matrix (e.g., normalized confusion matrix)
#' @param cutoff Threshold for binarizing the matrix (default: 0.1)
#' @param resolution Resolution parameter for Louvain clustering (default: 1.0)
#' @param algorithm Clustering algorithm: "louvain" or "leiden" (if igraph supports it)
#'
#' @return Integer vector of cluster assignments
#'
#' @details
#' The function:
#' 1. Binarizes the adjacency matrix using the cutoff
#' 2. Creates a graph from the binary matrix
#' 3. Applies Louvain/Leiden clustering to identify groups of connected clusters
#'
#' @examples
#' # Create a sample adjacency matrix
#' adj <- matrix(c(0, 0.3, 0.05, 0.3, 0, 0.02, 0.05, 0.02, 0), nrow = 3)
#' rownames(adj) <- colnames(adj) <- c("A", "B", "C")
#' cluster_adjacency_matrix(adj, cutoff = 0.1)
#'
#' @export
cluster_adjacency_matrix <- function(adj_matrix,
                                     cutoff = 0.1,
                                     resolution = 1.0,
                                     algorithm = "louvain") {

    adj_matrix <- as.matrix(adj_matrix)

    # Binarize the matrix
    binary_mat <- (adj_matrix > cutoff) * 1L

    # Ensure symmetric
    binary_mat <- pmax(binary_mat, t(binary_mat))
    diag(binary_mat) <- 0

    # Create graph
    g <- igraph::graph_from_adjacency_matrix(
        binary_mat,
        mode = "undirected",
        weighted = NULL,
        diag = FALSE
    )

    # Apply clustering
    if (algorithm == "louvain") {
        partition <- igraph::cluster_louvain(g, resolution = resolution)
    } else if (algorithm == "leiden") {
        if (requireNamespace("leiden", quietly = TRUE)) {
            partition <- leiden::leiden(binary_mat, resolution_parameter = resolution)
            # leiden package returns numeric vector directly
            return(as.integer(partition))
        } else {
            warning("leiden package not installed, using louvain instead")
            partition <- igraph::cluster_louvain(g, resolution = resolution)
        }
    } else {
        stop("algorithm must be 'louvain' or 'leiden'")
    }

    # Get membership
    groups <- igraph::membership(partition)

    return(as.integer(groups))
}

#' @title Merge clusters
#' @description Merge clusters based on group assignments
#'
#' @param labels Original cluster labels
#' @param groups New group assignments (from cluster_adjacency_matrix)
#' @param label_map Optional: named vector mapping old labels to new labels
#'
#' @return Factor with merged cluster labels
#'
#' @export
merge_clusters <- function(labels, groups, label_map = NULL) {
    labels <- as.character(labels)
    unique_labels <- sort(unique(labels))

    if (length(groups) != length(unique_labels)) {
        stop(sprintf("groups length (%d) doesn't match number of unique labels (%d)",
                     length(groups), length(unique_labels)))
    }

    # Create mapping from old labels to new group numbers
    if (is.null(label_map)) {
        label_map <- setNames(as.character(groups), unique_labels)
    }

    # Apply mapping
    new_labels <- label_map[labels]

    return(factor(new_labels))
}

#' @title Get connection matrix between clusterings
#' @description Compute connection matrix showing overlap between two clusterings
#'
#' @param labels1 First clustering labels (e.g., low resolution)
#' @param labels2 Second clustering labels (e.g., high resolution)
#' @param min_percent Minimum percentage threshold to consider a connection (default: 0.1)
#'
#' @return Binary connection matrix indicating which clusters from labels2 are connected
#'
#' @export
get_connection_matrix <- function(labels1, labels2, min_percent = 0.1) {
    # Count co-occurrences
    df <- data.frame(
        key1 = as.character(labels1),
        key2 = as.character(labels2),
        stringsAsFactors = FALSE
    )

    # Calculate counts and percentages
    count_df <- as.data.frame(table(df$key1, df$key2))
    names(count_df) <- c("key1", "key2", "counts")

    # Get cluster sizes for key2
    sizes <- table(df$key2)
    count_df$size <- sizes[as.character(count_df$key2)]
    count_df$percent <- count_df$counts / count_df$size

    # Filter by minimum percentage
    count_df <- count_df[count_df$percent > min_percent, ]

    # Find key1 groups that connect multiple key2 clusters
    key1_counts <- table(count_df$key1)
    connecting_key1 <- names(key1_counts)[key1_counts > 1]
    count_df <- count_df[count_df$key1 %in% connecting_key1, ]

    # Build connection matrix for key2 clusters
    unique_key2 <- sort(unique(labels2))
    n <- length(unique_key2)
    conn_mat <- matrix(0L, nrow = n, ncol = n)
    rownames(conn_mat) <- colnames(conn_mat) <- unique_key2

    # For each key1 that connects multiple key2 clusters, mark connections
    for (k1 in unique(count_df$key1)) {
        connected_k2 <- count_df$key2[count_df$key1 == k1]
        if (length(connected_k2) > 1) {
            # Create all pairs
            pairs <- utils::combn(as.character(connected_k2), 2)
            for (j in seq_len(ncol(pairs))) {
                i1 <- pairs[1, j]
                i2 <- pairs[2, j]
                conn_mat[i1, i2] <- conn_mat[i2, i1] <- 1L
            }
        }
    }

    return(conn_mat)
}

#' @title Bhattacharyya distance
#' @description Calculate Bhattacharyya distance between two probability distributions
#'
#' @param p First probability distribution
#' @param q Second probability distribution
#'
#' @return Bhattacharyya distance (non-negative value)
#'
#' @export
bhattacharyya_distance <- function(p, q) {
    # Handle zero-sum vectors
    p_sum <- sum(p)
    q_sum <- sum(q)
    
    if (p_sum == 0 || q_sum == 0) {
        warning("One of the distributions has zero sum")
        return(Inf)
    }
    
    # Ensure non-negative and normalize to sum to 1
    p <- pmax(p, 0)
    q <- pmax(q, 0)

    p <- p / p_sum
    q <- q / q_sum

    # Calculate Bhattacharyya coefficient: BC = sum(sqrt(p_i * q_i))
    # BC = 1 when distributions are identical
    # BC = 0 when distributions have no overlap
    bc <- sum(sqrt(p * q))

    # Calculate distance: D = -log(BC)
    # D = 0 when BC = 1 (identical)
    # D = Inf when BC = 0 (no overlap)
    dist <- -log(bc)

    # Handle NaN (shouldn't happen with proper input)
    if (is.nan(dist)) {
        return(0)
    }

    return(dist)
}

#' @title Bhattacharyya distance matrix
#' @description Compute pairwise Bhattacharyya distances between distributions
#'
#' @param prob_matrix Matrix where each column is a probability distribution
#' @param flags Optional logical vector to filter rows
#'
#' @return Distance matrix
#'
#' @export
bhattacharyya_matrix <- function(prob_matrix, flags = NULL) {
    n <- ncol(prob_matrix)
    dist_mat <- matrix(0, nrow = n, ncol = n)

    if (!is.null(colnames(prob_matrix))) {
        rownames(dist_mat) <- colnames(dist_mat) <- colnames(prob_matrix)
    }

    for (i in seq_len(n - 1)) {
        for (j in seq(i + 1, n)) {
            if (!is.null(flags)) {
                p <- prob_matrix[flags, i]
                q <- prob_matrix[flags, j]
            } else {
                p <- prob_matrix[, i]
                q <- prob_matrix[, j]
            }

            d <- bhattacharyya_distance(p, q)
            dist_mat[i, j] <- dist_mat[j, i] <- d
        }
    }

    return(dist_mat)
}

#' @title Get distance matrix between clusters
#' @description Compute pairwise distances between cluster centroids
#'
#' @param X Expression matrix (cells x features)
#' @param clusters Cluster assignments
#' @param labels Optional: specific labels to include (in order)
#' @param method Distance method: "euclidean", "manhattan", "cosine"
#'
#' @return Distance matrix between cluster centroids
#'
#' @export
get_distance_matrix <- function(X, clusters, labels = NULL, method = "euclidean") {
    clusters <- as.factor(clusters)

    if (is.null(labels)) {
        labels <- levels(clusters)
    }

    # Compute centroids
    centroids <- matrix(0, nrow = length(labels), ncol = ncol(X))
    rownames(centroids) <- labels

    for (i in seq_along(labels)) {
        cl <- labels[i]
        idx <- which(clusters == cl)
        if (length(idx) > 0) {
            if (inherits(X, "Matrix")) {
                centroids[i, ] <- Matrix::colMeans(X[idx, , drop = FALSE])
            } else {
                centroids[i, ] <- colMeans(X[idx, , drop = FALSE])
            }
        }
    }

    # Compute pairwise distances
    if (method == "cosine") {
        # Cosine distance = 1 - cosine similarity
        norm_centroids <- centroids / sqrt(rowSums(centroids^2))
        sim_mat <- tcrossprod(norm_centroids)
        dist_mat <- 1 - sim_mat
    } else {
        dist_mat <- as.matrix(stats::dist(centroids, method = method))
    }

    return(dist_mat)
}

#' @title Find high resolution clustering
#' @description Find a clustering resolution that ensures minimum cluster size
#'
#' @param object Seurat object
#' @param min_cells Minimum cells per cluster
#' @param max_resolution Maximum resolution to try
#' @param algorithm Clustering algorithm ("louvain" or "leiden")
#'
#' @return The Seurat object with appropriate resolution clustering
#'
#' @keywords internal
.find_high_resolution <- function(object,
                                  min_cells = 100,
                                  max_resolution = 4,
                                  algorithm = "louvain") {
    resolution <- max_resolution

    while (resolution > 0.5) {
        message(sprintf("Clustering with resolution: %.1f", resolution))

        # Run clustering
        if (algorithm == "leiden") {
            object <- Seurat::FindClusters(
                object,
                resolution = resolution,
                algorithm = 4,  # Leiden
                verbose = FALSE
            )
        } else {
            object <- Seurat::FindClusters(
                object,
                resolution = resolution,
                algorithm = 1,  # Louvain
                verbose = FALSE
            )
        }

        # Store with resolution-specific name
        key <- sprintf("clustering_res%.1f", resolution)
        object[[key]] <- Seurat::Idents(object)

        # Check minimum cluster size
        cluster_sizes <- table(Seurat::Idents(object))
        if (min(cluster_sizes) >= min_cells) {
            break
        }

        resolution <- resolution - 0.5
    }

    return(object)
}
