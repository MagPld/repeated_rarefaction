test_threshold <- function(physeq, repeats = 10, t_min = 5, t_max = 250, t_step = 1, method = "NMDS", sample_id, groupb) {
  thresholds <- as.integer(seq.int(t_min, t_max, by = t_step))
  index_data <- matrix(nrow = 0, ncol = 3)
  colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")

  # Loops through all thresholds set in the parametteser above
  for (y in repeats) {
    for (x in thresholds) {
      print(paste("Running with ", y, " repeats and ", x, " threshold"))
      index <- get_index(physeq, sample_data(physeq), y, threshold = x, method, sample_id, groupb)
      print("test5")
      index_data <- rbind(index_data, c(toString(y), x, index))
      print("test2")
    }
  }


  # colnames(index_data) <- c("Repeat_Amount", "Threshold", "Index")
  index_data <- as.data.frame(index_data)
  index_data$Threshold <- as.numeric(index_data$Threshold)
  index_data$Index <- as.numeric(index_data$Index)

  ggplot(index_data, aes(x = Threshold, y = Index, color = Repeat_Amount)) +
    geom_point() +
    labs(x = "Rarefaction Threshold", y = "Calinski-Harabasz pseudo F-statistic") +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F)
}

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


  clusterSim::index.G1(just_positions, cluster_true) #Throws error about matrix (invalid nrow value)
  print("test4")
  return(clusterSim::index.G1(just_positions, cluster_true)) #Throws error: can't find index.G1 function
}
