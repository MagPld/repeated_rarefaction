#' This function is a wrapper designed around "repeated_rarefaction". It allows.
#'
#' @param input A phyloseq object or a table of OTU/ASV counts.
#' @param repeats A positive integer. Indicates the amount of repeats run. A value = 1 means no repeats are used.
#' @param t_min An integer. Minimum value for the threshold testing range.
#' @param t_max An integer. Maximum value for the threshold testing range.
#' @param t_step An float. Step value for the threshold testing range. A value between 0 and 1 will cause the same threshold to be tested more than once.
#' @param sample_id A string. Name of the column with sample IDs.
#' @param group A string. Name of the column with data to group the points by.
#' @param cores An int. The amount of cores used when running. The default = 4.
#' @returns A graph with rarefaction threshold as the x axis, and the y axis is the calculated Calinski-Harabasz F-statistic. Each dot indicates the clustering score for the repeated ordination on a single threshold value. The higher score the better.
#' @examples
#' test_threshold(HLCYG_physeq_data, repeats = 10, t_min = 200, t_max =1500, t_step = 5, group="location")
test_threshold <- function(input, repeats = 10, t_min = 50, t_max = 250, t_step = 1, group = "sample_id", cores = 4) {
  # Check if input is a Phyloseq object

  if (inherits(input, "phyloseq")) {
    physeq <- input
  } else {
    stop("Input must be a Phyloseq object, including a count table and sample data.")
  }

  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable.
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))
  
  # Determine the thresholds to loop over between min and max and specified step.
  thresholds <- seq.int(t_min, t_max, by = t_step)
  # Create the plot list which will hold the ordinations and the matrix which 
  # will hold index calculations results.
  plots <- list()
  index_data <- matrix(nrow = 0, ncol = 3)
  
  for (y in repeats) {
    # Create lists to store ordination points coordinates for a specific repeat 
    # amount and centroid coordinates.
    ordination_points <- list()
    consensus_coordinates <- list()
    
    for (x in thresholds) {
      message(paste("Running with", y, "repeats and", x, "threshold"))

      step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold = x, repeats = y)
      step2 <- ord_and_mean(step1$rarefied_matrix_list, repeats)

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
    
    # Flatten the list of lists containing individual repeat attempts for each threshold
    # After flattening, each element in flattened_list corresponds to a single threshold and
    # is a combined data frame (or matrix) of coordinates from all repeats for that threshold.
    flattened_list <- lapply(ordination_points, function(sublist) {
      do.call(rbind, sublist)
    })
    
    # =========================
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
      # 'consensus_coordinates' should also be normalized.
      norm_consensus <- normalize_coords(consensus_coordinates[[current_key]], axis_means, axis_sds)

      step3 <- plot_rep_raref(split_reps, norm_consensus, sample_data(physeq), color = group, group = group, cloud = TRUE, ellipse = TRUE, title = paste0("Threshold: ", x))
      plot <- step3$plot

      plots[[paste0("repeat_number ", y)]][[current_key]] <- plot + xlim(x_limits) + ylim(y_limits)
      
      # ======================== Index calculation
      info <- sample_data(physeq)

      # Fix the normalized consensus coordinates adding group variable
      norm_consensus_df <- data.frame
      temp_df <- as.data.frame(norm_consensus)
      temp_df$sample_id <- rownames(split_reps[[1]])
      rownames(temp_df) <- rownames(split_reps[[1]])
      temp_df[[group]] <- info[[group]]
      norm_consensus_df <- rbind(temp_df)

      # Extract just position data for the individual points
      just_positions_consensus <- norm_consensus_df[,1:2]

      # Determine clusters according to the "group" variable
      clusters_consensus <- unlist(norm_consensus_df[[group]])
      existing_groups_consensus <- unique(clusters_consensus)

      # Convert cluster labels to numbers
      cluster_consensus_numbers <- match(clusters_consensus, existing_groups_consensus)

      # Calculate index
      index <- clusterSim::index.G1(just_positions_consensus, cluster_consensus_numbers)

      index_data <- rbind(index_data, c(toString(y), x, index))
    }
    
  }
  
  # Format the index_data matrix for plotting
  index_data <- as.data.frame(index_data)
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)
  # Plot the indexes
  index_plot <- ggplot(index_data, aes(x = Threshold, y = Index, color = Repeat_Amount)) +
    geom_point() +
    labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
    geom_smooth(method = "loess", formula = 'y ~ x')
  
  output <- list(index_plot = index_plot, ordination_plots = plots)
  class(output) <- "test_threshold"
  return(output)
}

print.test_threshold <- function(x, ...) {
  # Print only the index plot
  print(x[["index_plot"]])
  invisible(x)  # standard practice for print methods
}
