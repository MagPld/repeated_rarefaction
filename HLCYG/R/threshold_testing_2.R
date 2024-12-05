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
  
  # a group information is not normally necessary for rep_raref but it is when is run to calculate
  # the index for clustering performance
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
  
  # The only thing I need theoretically is true position and cluster membership. 
  # However this still misses something. Because since we are not calculating cluster dispersion
  # or not taking it into account, then what happens is that the median positions cannot
  # worsen that much
  
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
  # It's possible to change which index is used
  # Currently it's set to G1
  cluster_numbers <- as.numeric(unlist(cluster_numbers))

  return(clusterSim::index.G1(just_positions, cluster_numbers))
}
