#' Test a range of threshold values for repeated rarefaction.
#'
#' @param physeq A phyloseq object.
#' @param repeats A positive integer. Indicates the amount of repeats run. A value = 1 means no repeats are used.
#' @param t_min An integer. Minimum value for the threshold testing range.
#' @param t_max An integer. Maximum value for the threshold testing range.
#' @param t_step An float. Step value for the threshold testing range. A value between 0 and 1 will cause the same threshold to be tested more than once.
#' @param method A string. Currently only "NMDS" supported.
#' @param sample_id A string. Name of the column with sample IDs.
#' @param groupb A string. Name of the column with data to group the points by.
#' @returns A graph with rarefaction threshold as the x axis, and the y axis is the calculated Calinski-Harabasz F-statistic. Each dot is the score for one threshold attempt. The higher score the better.
#' @examples
#' test_threshold(HLCYG_physeq_data, repeats = 5, t_step = 10, method= "NMDS", sample_id = "sample_id", groupb = "location")
test_threshold_2<- function(physeq, repeats = 10, t_min = 50, t_max = 250, t_step = 1, group="sample_id") {
  
  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable.
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))

  if (!(is.double(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }

  # This introduces an error if I pass a vector of repeats instead than a single element.
  # if (!(repeats == round(repeats))){
  #   stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  # }

  thresholds <- as.integer(seq.int(t_min, t_max, by = t_step))
  index_data <- matrix(nrow = 0, ncol = 3)
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")

  # Loops through all thresholds set in the parametteser above
  for (y in repeats) {
    for (x in thresholds) {
      print(paste("Running with ", y, " repeats and ", x, " threshold"))
      index <- get_index_2(physeq, sample_data(physeq), repeats= y, threshold = x, group = group)
      index_data <- rbind(index_data, c(toString(y), x, index))
    }
  }
  
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data <- as.data.frame(index_data)
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)

  ggplot(index_data, aes(x = Threshold, y = Index, color = Repeat_Amount)) +
    geom_point() +
    labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F)
}


#' Performs repeatead rarefaction once and calculates a performance value
#'
#' @param physeq A phyloseq object.
#' @param sample_info Sample Data object from the phyloseq package. Contains data info without repetitions.
#' @param repeats A positive integer. Indicates the amount of repeats run. A value = 1 means no repeats are used.
#' @param threshold A positive integer. The threshold value for the rarefaction.
#' @param method A string. Currently only "NMDS" supported.
#' @param sample_id A string. Name of the column with sample IDs.
#' @param groupb A string. Name of the column with data to group the points by.
#' @returns Returns the Calinski-Harabasz F-statistic score.
get_index_2 <- function(physeq, sample_info, repeats, threshold, group) {
  # Calls the first two steps of the repeated rarefaction methods

  step1 <- rep_raref_2(data.frame(t(otu_table(physeq))), threshold, repeats)
  step2 <- ord_and_mean_2(step1$rarefied_matrix_list, repeats)
  
  # Combine aligned ordinations into one data frame for calculation
  aligned_df <- data.frame()
  
  # Setup needed resources
  aligned_ordinations <- step2$aligned_ordinations
  info <- sample_data(physeq)
  
  for (i in 1:length(aligned_ordinations)) {
    temp_df <- as.data.frame(aligned_ordinations[[i]])
    colnames(temp_df) <- c("Dim1", "Dim2")
    temp_df$sample_id <- rownames(aligned_ordinations[[i]])
    temp_df$ordination <- paste0("Ordination", i)
    temp_df[[group]] <- info[[group]]
    aligned_df <- rbind(aligned_df, temp_df)
  }

  # Extract just position data
  just_positions <- aligned_df[1:2]
  
  # Determine clusters according to the "group" variable
  clusters <- unlist(aligned_df[[group]])
  # Identify unique values (corresponsing to number of clusters)
  existing_groups <- unique(clusters)
  
  # Index function requires clusters to be numbered.
  # Therefore whatever the cluster names are called, they are here converted to numbers
  cluster_numbers <- list()
  for (x in 1:length(clusters)) {
    cluster_numbers[x] <- which(existing_groups == clusters[x])
  }

  # Calculate index
  cluster_numbers <- as.numeric(unlist(cluster_numbers))
  ## Auto calculation of centroids
  #index <- clusterSim::index.DB(just_positions, cluster_numbers)
  
  ## Manual caculation with previously mshape determined centroids

  consensus_coords <- step2$consensus_coordinates
  
  # Calculate intra-cluster spread
  intra_cluster_spread <- sapply(unique(cluster_numbers), function(k) {
    cluster_points <- just_positions[cluster_numbers == k, ]  # Points in cluster k
    centroid <- consensus_coords[k, ]  # Precomputed centroid for cluster k
    mean(rowSums((cluster_points - centroid)^2))  # Mean squared Euclidean distance
  })
  
  # Calculate inter-cluster distances
  inter_cluster_distances <- as.matrix(dist(consensus_coords))  # Pairwise distances between centroids
  diag(inter_cluster_distances) <- Inf  # Ignore self-distances
  
  # Calculate Davies-Bouldin Index
  db_index <- mean(sapply(unique(cluster_numbers), function(i) {
    s_i <- intra_cluster_spread[i]
    ratios <- sapply(setdiff(unique(cluster_numbers), i), function(j) {
      s_j <- intra_cluster_spread[j]
      m_ij <- inter_cluster_distances[i, j]
      (s_i + s_j) / m_ij  # Ratio for cluster i and j
    })
    max(ratios)  # Take the maximum ratio for cluster i
  }))

  return(db_index)
}

test_threshold_test <- function(input, repeats = 10, t_min = 50, t_max = 250, t_step = 1, group = "sample_id") {
  # Check if input is a Phyloseq object or a file path
  if (inherits(input, "phyloseq")) {
    physeq <- input
  } else if (is.character(input) && file.exists(input)) {
    count_table <- read.table(input, header = TRUE, row.names = 1, check.names = FALSE)
    physeq <- phyloseq(otu_table(count_table, taxa_are_rows = TRUE))
  } else {
    stop("Input must be a Phyloseq object or a file path to an OTU count table (.csv or .tsv).")
  }
  
  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable.
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))
  
  thresholds <- seq.int(t_min, t_max, by = t_step)
  
  # Store ordination points and centroids
  ordination_points <- list()
  consensus_coordinates <- list()
  
  for (y in repeats) {
    for (x in thresholds) {
      message(paste("Running with", y, "repeats and", x, "threshold"))
      
      step1 <- rep_raref_2(data.frame(t(otu_table(physeq))), threshold = x, repeats = y)
      step2 <- ord_and_mean_2(step1$rarefied_matrix_list, repeats)
      
      # step2 should return something like: 
      # $aligned_ordinations: A list of ordination coordinates for each repeat
      # $consensus_coordinates: a matrix or df of centroid coords
      
      # Store aligned ordination points for this threshold
      # Note: They may not actually be "aligned" in any external sense now, 
      # but are just coordinates from each repetition's ordination.
      ordination_points[[paste0("threshold_", x)]] <- step2$aligned_ordinations
      
      # Store consensus coordinates for this threshold
      consensus_coords <- step2$consensus_coordinates
      colnames(consensus_coords) <- c("Dim1", "Dim2")
      consensus_coordinates[[paste0("threshold_", x)]] <- consensus_coords
    }
  }
  
  # Flatten the list of lists containing individual repeat attempts for each threshold
  # After flattening, each element in flattened_list corresponds to a single threshold and
  # is a combined data frame (or matrix) of coordinates from all repeats for that threshold.
  flattened_list <- lapply(ordination_points, function(sublist) {
    do.call(rbind, sublist)
  })
  
  # =========================
  # Remove Procrustes alignment
  # Instead, we will normalize all coordinates to unit variance.
  
  # Combine all points across all thresholds to compute global scaling parameters
  all_points <- do.call(rbind, flattened_list)
  
  # Compute global means and sds for each axis
  axis_means <- colMeans(all_points, na.rm = TRUE)
  axis_sds   <- apply(all_points, 2, sd, na.rm = TRUE)
  
  # Define a normalization function
  normalize_coords <- function(mat, means, sds) {
    sweep(sweep(mat, 2, means, "-"), 2, sds, "/")
  }
  
  # Normalize each threshold's ordination points
  scaled_points <- lapply(flattened_list, function(x) normalize_coords(x, axis_means, axis_sds))
  
  # Determine global axis limits after normalization (with some tolerance)
  all_scaled_points <- do.call(rbind, scaled_points)
  x_limits <- range(all_scaled_points[, 1], na.rm = TRUE)
  y_limits <- range(all_scaled_points[, 2], na.rm = TRUE)
  x_limits <- c(x_limits[1] - 0.2, x_limits[2] + 0.2)
  y_limits <- c(y_limits[1] - 0.2, y_limits[2] + 0.2)
  
  # =====================
  # Plotting and Index Calculation
  
  plots <- list()
  index_data <- matrix(nrow = 0, ncol = 4)
  
  for (y in repeats) {
    for (x in thresholds) {
      
      # Retrieve normalized data for this threshold
      # scaled_points[[paste0("threshold_", x)]] is a combined matrix with all repeats
      # We need to split by repeats again
      current_key <- paste0("threshold_", x)
      df <- as.data.frame(scaled_points[[current_key]])
      
      # Split into repeats
      n_points <- nrow(df)
      # each threshold run had 'repeats' ordinations, so we split equally
      split_reps <- split(df, ceiling(seq_len(n_points) / (n_points/y)))
      
      # Plot
      # Now that we aren't using procrustes, 'consensus_coordinates' should also be normalized.
      # Let's normalize them so they match the scaling of points.
      # NOTE: We must apply the same normalization to the consensus coords:
      norm_consensus <- normalize_coords(consensus_coordinates[[current_key]], axis_means, axis_sds)
      
      step3 <- plot_rep_raref_2(split_reps, norm_consensus, sample_data(physeq), color = group, group = group, cloud = TRUE, ellipse = TRUE)
      plot <- step3$plot
      plots[[current_key]] <- plot + xlim(x_limits) + ylim(y_limits)
      
      # ======================== Index calculation
      
      # aligned_df <- data.frame()
      info <- sample_data(physeq)
      
      # Fix the repetitions points and df
      # for (i in 1:length(split_reps)) {
      #   temp_df <- as.data.frame(split_reps[[i]])
      #   colnames(temp_df) <- c("Dim1", "Dim2")
      #   temp_df$sample_id <- rownames(split_reps[[i]])
      #   temp_df$ordination <- paste0("Ordination", i)
      #   temp_df[[group]] <- info[[group]]
      #   aligned_df <- rbind(aligned_df, temp_df)
      # }
      
      # Fix the normalized consensus coordinates adding group variable
      norm_consensus_df <- data.frame
      temp_df <- as.data.frame(norm_consensus)
      temp_df$sample_id <- rownames(split_reps[[1]])
      rownames(temp_df) <- rownames(split_reps[[1]])
      temp_df[[group]] <- info[[group]]
      norm_consensus_df <- rbind(temp_df)
      
      # Extract just position data for the individual points
      # just_positions <- aligned_df[, c("Dim1", "Dim2")]
      just_positions_consensus <- norm_consensus_df[,1:2]
      
      # Determine clusters according to the "group" variable (for the points)
      # clusters <- unlist(aligned_df[[group]])
      # existing_groups <- unique(clusters)
      
      # Determine clusters for the consensus coordinates
      clusters_consensus <- unlist(norm_consensus_df[[group]])
      existing_groups_consensus <- unique(clusters_consensus)
      
      # Convert cluster labels to numbers
      # cluster_numbers <- match(clusters, existing_groups)
      cluster_consensus_numbers <- match(clusters_consensus, existing_groups_consensus)
      
      # Calculate index
      index <- clusterSim::index.G1(just_positions_consensus, cluster_consensus_numbers)
      # index_consensus <- clusterSim::index.G1(just_positions_consensus, cluster_consensus_numbers)
      
      index_data <- rbind(index_data, c(toString(y), x, index))

    }
  }
  
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data <- as.data.frame(index_data)
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)
  # index_data$Index_consensus <- as.numeric(index_data$Index_consensus)
  
  # Plot the indexes
  index_plot <- ggplot(index_data, aes(x = Threshold, y = Index, color = Repeat_Amount)) +
    geom_point() +
    labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)
  
  # index_consensus_plot <- ggplot(index_data, aes(x = Threshold, y = Index_consensus, color = Repeat_Amount)) +
  #   geom_point() +
  #   labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
  #   geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)
  
  return(list(index_plot = index_plot, ordination_plots = plots))
}
