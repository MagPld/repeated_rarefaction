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
test_threshold <- function(physeq, repeats = 10, t_min = 5, t_max = 250, t_step = 1, method = "NMDS", sample_id, groupb) {
  thresholds <- as.integer(seq.int(t_min, t_max, by = t_step))
  index_data <- matrix(nrow = 0, ncol = 3)
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")

  # Loops through all thresholds set in the parametteser above
  for (y in repeats) {
    for (x in thresholds) {
      print(paste("Running with ", y, " repeats and ", x, " threshold"))
      index <- get_index(physeq, sample_data(physeq), repeats= y, threshold = x, method, sample_id, groupb)
      index_data <- rbind(index_data, c(toString(y), x, index))
    }
  }


  # colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data <- as.data.frame(index_data)
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)

  print(index_data)

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
get_index <- function(physeq, sample_info, repeats, threshold, method, sample_id, groupb) {
  # Calls the first two steps of the repeated rarefaction methods
  step1 <- rep_raref(as(otu_table(physeq), "matrix"), sample_data(physeq), threshold, repeats)
  step2 <- ord_and_mean(step1$repeat_count, step1$repeat_info, sample_info, repeats, method, sample_id)

  # Extract just position data and create a list of "true" clusters
  just_positions <- as.data.frame(step2$df_median[,2:3])
  cluster_true_pre <- unlist(step2$df_median[groupb])
  cluster_names <- list()
  cluster_true <- cluster_true_pre

  # Index function requires clusters to be numbered.
  # Therefore whatever the cluster names are called, they are here converted to numbers
  existing_groups <- unique(cluster_true)
  for (x in 1:length(cluster_true_pre)) {
    cluster_true[x] <- which(existing_groups == cluster_true_pre[x])
  }

  # Calculate index
  # It's possible to change which index is used
  # Currently it's set to G1
  cluster_true <- as.numeric(unlist(cluster_true))

  return(clusterSim::index.G1(just_positions, cluster_true))
}
