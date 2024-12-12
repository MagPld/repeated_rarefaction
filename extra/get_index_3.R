get_index_3 <- function(physeq, sample_info, repeats, threshold, group) {
  # Calls the first two steps of the repeated rarefaction methods
  
  # Step 1: Perform rarefaction and ordination
  step1 <- rep_raref_2(data.frame(t(otu_table(physeq))), threshold, repeats)
  step2 <- ord_and_mean_2(step1$rarefied_matrix_list, repeats)
  
  # Combine aligned ordinations into one data frame
  info <- sample_data(physeq)
  aligned_df <- do.call(rbind, lapply(1:length(step2$aligned_ordinations), function(i) {
    temp_df <- as.data.frame(step2$aligned_ordinations[[i]])
    colnames(temp_df) <- c("Dim1", "Dim2")
    temp_df$sample_id <- rownames(temp_df)
    temp_df$ordination <- paste0("Ordination", i)
    temp_df[[group]] <- info[[group]]
    temp_df
  }))
  
  # Extract position data as a matrix for efficiency
  just_positions <- as.matrix(aligned_df[, c("Dim1", "Dim2")])
  
  # Map unique cluster names to numeric labels
  clusters <- aligned_df[[group]]
  cluster_numbers <- match(clusters, unique(clusters))
  
  # Use precomputed centroids from step2
  consensus_coords <- as.matrix(step2$consensus_coordinates)
  
  # Calculate intra-cluster spread
  intra_cluster_spread <- sapply(unique(cluster_numbers), function(k) {
    cluster_points <- just_positions[cluster_numbers == k, ]
    centroid <- consensus_coords[k, ]
    mean(rowSums((cluster_points - centroid)^2))  # Mean squared Euclidean distance
  })
  
  # Precompute inter-cluster distances
  inter_cluster_distances <- as.matrix(dist(consensus_coords))
  diag(inter_cluster_distances) <- Inf  # Ignore self-distances
  
  # Parallelize Davies-Bouldin Index calculation
  cl <- makeCluster(detectCores() - 1)  # Use all but one core
  registerDoParallel(cl)
  
  db_index <- mean(foreach(i = unique(cluster_numbers), .combine = c, .packages = "base") %dopar% {
    s_i <- intra_cluster_spread[i]
    ratios <- sapply(setdiff(unique(cluster_numbers), i), function(j) {
      s_j <- intra_cluster_spread[j]
      m_ij <- inter_cluster_distances[i, j]
      (s_i + s_j) / m_ij
    })
    max(ratios)  # Take the maximum ratio for cluster i
  })
  
  stopCluster(cl)  # Shut down the cluster
  
  return(db_index)
}