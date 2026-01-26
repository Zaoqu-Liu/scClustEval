# Tests for clustering functions

test_that("cluster_adjacency_matrix works correctly", {
    skip_if_not_installed("igraph")

    # Create adjacency matrix with clear groups
    adj <- matrix(0, nrow = 4, ncol = 4)
    rownames(adj) <- colnames(adj) <- c("A", "B", "C", "D")

    # Group 1: A-B connected
    adj["A", "B"] <- adj["B", "A"] <- 0.5

    # Group 2: C-D connected
    adj["C", "D"] <- adj["D", "C"] <- 0.5

    groups <- cluster_adjacency_matrix(adj, cutoff = 0.1)

    expect_equal(length(groups), 4)
    # A and B should be in same group
    expect_equal(groups[1], groups[2])
    # C and D should be in same group
    expect_equal(groups[3], groups[4])
    # Groups should be different
    expect_true(groups[1] != groups[3])
})

test_that("cluster_adjacency_matrix handles isolated nodes", {
    skip_if_not_installed("igraph")

    adj <- matrix(0, nrow = 3, ncol = 3)
    rownames(adj) <- colnames(adj) <- c("A", "B", "C")

    # Only A-B connected
    adj["A", "B"] <- adj["B", "A"] <- 0.5

    groups <- cluster_adjacency_matrix(adj, cutoff = 0.1)

    expect_equal(length(groups), 3)
    expect_equal(groups[1], groups[2])  # A and B together
    expect_true(groups[3] != groups[1])  # C isolated
})

test_that("get_connection_matrix identifies correct connections", {
    # Low resolution: 2 clusters
    labels1 <- c(rep("L1", 50), rep("L2", 50))

    # High resolution: 4 clusters (2 sub-clusters per low-res cluster)
    labels2 <- c(rep("H1", 25), rep("H2", 25), rep("H3", 25), rep("H4", 25))

    conn_mat <- get_connection_matrix(labels1, labels2)

    expect_equal(dim(conn_mat), c(4, 4))

    # H1 and H2 should be connected (both from L1)
    expect_equal(conn_mat["H1", "H2"], 1)

    # H3 and H4 should be connected (both from L2)
    expect_equal(conn_mat["H3", "H4"], 1)

    # H1 and H3 should NOT be connected (different L clusters)
    expect_equal(conn_mat["H1", "H3"], 0)
})

test_that("merge_clusters works correctly", {
    labels <- factor(c("A", "A", "B", "B", "C", "C"))
    groups <- c(1, 1, 2)  # A and B merge into group 1, C stays as 2

    new_labels <- merge_clusters(labels, groups)

    expect_equal(length(new_labels), 6)
    # First 4 cells should be in same group
    expect_true(length(unique(new_labels[1:4])) == 1)
})

test_that("bhattacharyya_distance computes correctly", {
    # Identical distributions should have distance 0
    p <- c(0.25, 0.25, 0.25, 0.25)
    q <- c(0.25, 0.25, 0.25, 0.25)

    dist <- bhattacharyya_distance(p, q)
    expect_equal(dist, 0, tolerance = 1e-10)

    # Similar distributions should have small distance
    p <- c(0.4, 0.3, 0.2, 0.1)
    q <- c(0.35, 0.35, 0.2, 0.1)

    dist <- bhattacharyya_distance(p, q)
    expect_true(dist >= 0)
    expect_true(dist < 1)  # Similar distributions have small distance

    # Different distributions should have larger distance
    p <- c(0.9, 0.05, 0.03, 0.02)
    q <- c(0.02, 0.03, 0.05, 0.9)

    dist2 <- bhattacharyya_distance(p, q)
    expect_true(dist2 > dist)  # More different = larger distance
})

test_that("bhattacharyya_matrix is symmetric", {
    prob_mat <- matrix(c(
        0.7, 0.2, 0.1,
        0.1, 0.8, 0.1,
        0.2, 0.1, 0.7
    ), nrow = 3, byrow = TRUE)

    dist_mat <- bhattacharyya_matrix(t(prob_mat))

    # Should be symmetric
    expect_equal(dist_mat, t(dist_mat))

    # Diagonal should be 0
    expect_equal(diag(dist_mat), c(0, 0, 0))
})

test_that("get_distance_matrix computes euclidean correctly", {
    # Simple 2D case
    X <- matrix(c(
        0, 0,  # Cell 1, Cluster A
        1, 0,  # Cell 2, Cluster A
        10, 0,  # Cell 3, Cluster B
        11, 0   # Cell 4, Cluster B
    ), nrow = 4, byrow = TRUE)

    clusters <- factor(c("A", "A", "B", "B"))

    dist_mat <- get_distance_matrix(X, clusters)

    expect_equal(dim(dist_mat), c(2, 2))
    # Diagonal should be 0 (ignore names)
    expect_equal(as.numeric(diag(dist_mat)), c(0, 0))

    # Distance between cluster centroids
    # Centroid A = (0.5, 0), Centroid B = (10.5, 0)
    # Distance = 10
    expect_equal(dist_mat["A", "B"], 10, tolerance = 1e-10)
})
