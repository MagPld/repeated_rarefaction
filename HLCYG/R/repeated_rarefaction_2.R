#' Perform one iteration of repeated rarefaction and produce ordination plot
#'
#' @param physeq A phyloseq object.
#' @param repeats A positive integer. Indicates the amount of repeats run.
#' A value = 1 means no repeats are used.
#' @param threshold A positive integer. The threshold value for the rarefaction.
#' @param method A string. Currently only "NMDS" supported.
#' @param sample_id A string. Name of the column with sample IDs.
#' @param colorb A string. Name of the column with data to colour the graph by.
#' @param shapeb A string. Name of the column with data to shape the points on the graph by.
#' @param cloud Boolean. Aesthetic setting for the graph. Default is FALSE.
#' TRUE shows datapoints for all repeats.
#' @param ellipse Boolean. Aesthetic setting for the graph. Default is TRUE.
#' @returns A list with repeat count table, repeat info table, ordination object,
#' physeq object, dataframe with all repeat ordination positions, dataframe with median ordination positions, and ordination plot.
#' @examples
#' repeated_rarefaction(HLCYG_physeq_data, repeats=5, threshold=500, method="NMDS", colorb="sample_id", shapeb="location", T, F)
#' repeated_rarefaction(HLCYG_physeq_data, repeats=10, threshold=250, method="NMDS", colorb="sample_id", shapeb="location", T, T)
repeated_rarefaction_2 <- function(physeq, repeats = 50, threshold = 250, colorb="sample_id", shapeb="sample_id", cloud = FALSE, ellipse = TRUE) {
  
  # Make the rownames of the Phyloseq object a new "sample_id" variable for the sample data.
  # (this covers the case in which no sample_id column is present in the sample data)
  # Then set it to a separate variable.
  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))
  
  # ============ Checks and warnings
  
  if (!(colorb %in% names(sample_data(physeq)))) {
    stop(paste("'",colorb,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to color the ordination plot by.", sep=""))
  }
  if (!(shapeb %in% names(sample_data(physeq)))) {
    stop(paste("'",shapeb,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to shape points in the ordination plot by.", sep=""))
  }

  if (!(is.double(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }

  if (!(repeats == round(repeats))){
    stop(paste("Input for repeats: '", repeats, "' is not an integer.", sep=""))
  }

  if (!(is.double(threshold))){
    stop(paste("Input for threshold: '" ,threshold, "' is not an integer.", sep=""))
  }

  if (repeats <=4 & ellipse == TRUE){
    warning("Too few repeats to draw confidence ellipses.")
    ellipse <- F
  }
  
  # Perform the different steps of the repeated rarefaction algorithm 
  step1 <- rep_raref(data.frame(t(otu_table(physeq))), threshold, repeats)
  step2 <- ord_and_mean(step1$rarefied_matrix_list, repeats)
  step3 <- plot_rep_raref(step2$aligned_ordinations, step2$consensus_coordinates, colorb, shapeb, cloud, ellipse)
  print(step3)
  return(invisible(list("repeat_count" = step1$repeat_count, "repeat_info" = step1$repeat_info, "ordinate_object" = step2$ordinate_object, "physeq_object" = step2$physeq_object, "df_all" = step2$df_all, "df_median" = step2$df_median, "plot" = step3)))
}

#' Rarefaction is performed repeatedly depending on input and a matrix containing all data is created together with a an info file reflecting it.
#'
#' @param count A matrix. An otu-table containing the count data the rarefaction should be performed on.
#' @param info Sample Data object from the phyloseq package.
#' Should contain info about the count data with matching column names.
#' @param threshold An integer. The threshold value at which to perform rarefaction.
#' @param repeats An integer. The amount of repeated rarefactions to perform.
#' A value = 1 means only one iteration of rarefaction is perfomed and therefore no repeats.
#' @returns A list containing a matrix with the repeated count table and another
#' matrix with the repeated info.
rep_raref <- function(count, threshold, repeats) {
  if (repeats == 0) {
    warning("repeats can't be 0. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  if (repeats < 0) {
    warning("repeats can't be negative. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  # Set up working files
  rarefied_matrices <- list()

  # Perform repeated rarefaction and store the normalized results in a list
  # TODO: check if the use of bray curtis actually require normalization
  if (repeats >= 1) {
    for (i in 1:repeats) {
      rarefied_count <- rrarefy(count, sample = threshold)
      rarefied_matrices[[i]] <- rarefied_count
    }
  }

  return(list("rarefied_matrix_list"=rarefied_matrices))
}

#' Calculates ordination on the repeated count input and as well as a median result for each original datapoint.
#'
#' @param repeat_count A matrix. Contains the repeated count data.
#' @param repeat_info A matrix. Contains the repeated data info.
#' @param sample_info Sample Data object from the phyloseq package. Contains data info without repetitions.
#' @param repeats An integer. The amount of repeats that has been performed on the data.
#' @param method A string. Currently only "NMDS" supported.
#' @param sample_id A string. Name of the column with sample IDs.
#' @returns Returns a list containg an ordination object, a phyloseq object with
#' the repeated count data, a dataframe with all positions from the ordination
#' calculaction, and a datafram with just the median position from the calculation.
ord_and_mean <- function(rarefied_matrix_list, repeats) {
  
  #========================= ordinations and plots generation
  
  # Initialize a list to store ordinations
  ordinations <- list()
  
  # Set up parallel backend to use available cores
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Perform parallel computation using foreach
  results <- foreach(i = 1:length(rarefied_matrix_list), .packages = c('vegan', 'ggplot2')) %dopar% {
    # Calculate Bray-Curtis distance
    dist_matrix <- vegdist(rarefied_matrix_list[[i]], method = "bray")
    
    # Perform PCoA (Principal Coordinates Analysis)
    ordination <- cmdscale(dist_matrix, k = 2)
    
    # Convert ordination result to data frame
    ord_df <- as.data.frame(ordination)
    colnames(ord_df) <- c("PCoA1", "PCoA2")
    ord_df$Sample <- rownames(ordination)
    
    # Return a list containing the ordination and the plot
    list(ordination = ordination)
    
  }
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  
  # Extract ordinations and plots from the results
  for (i in 1:length(results)) {
    ordinations[[i]] <- results[[i]]$ordination
  }
  
  #================================ procrustes
  
  # Perform Procrustes analysis to align all ordinations to the first one
  aligned_ordinations <- lapply(ordinations, function(x) procrustes(ordinations[[1]], x)$Yrot)
  
  # Convert list to array for consensus calculation
  aligned_array <- array(unlist(aligned_ordinations), dim = c(nrow(aligned_ordinations[[1]]), ncol(aligned_ordinations[[1]]), length(aligned_ordinations)))
  
  # Compute consensus using mean shape
  consensus_coords <- mshape(aligned_array)

  return(list("aligned_ordinations" = aligned_ordinations, "consensus_coordinates" = consensus_coords))
}


#' Create a plot for the repeated rarefaction.
#'
#' @param ordination An ordination object.
#' @param physeq A phyloseq object containing the data to be plotted.
#' @param all_positions A dataframe containing position data for all repeats.
#' @param mediant_positions A dataframe containing median position data for each original data point.
#' @param color A string. Name of the column to color points by.
#' @param shape A string. Name of the column to shape points by.
#' @param cloud Boolean. Aesthetic setting for the graph. Default is FALSE.
#' TRUE shows datapoints for all repeats.
#' #' @param ellipse Boolean. Aesthetic setting for the graph. Default is TRUE.
#' @returns Returns an ordination plot.
plot_rep_raref <- function(aligned_ordinations, consensus_coordinates, color, shape, cloud, ellipse) {
  
  # Combine aligned ordinations into one data frame for plotting
  aligned_df <- data.frame()
  
  for (i in 1:length(aligned_ordinations)) {
    temp_df <- as.data.frame(aligned_ordinations[[i]])
    colnames(temp_df) <- c("Dim1", "Dim2")
    temp_df$Sample <- rownames(aligned_ordinations[[i]])
    temp_df$Ordination <- paste0("Ordination", i)
    aligned_df <- rbind(aligned_df, temp_df)
  }
  
  # Convert consensus_coordinates to a data frame
  consensus_df <- as.data.frame(consensus_coordinates)
  colnames(consensus_df) <- c("Dim1", "Dim2")
  consensus_df$Sample <- rownames(consensus_coordinates)
  
  # Plot all aligned ordinations
  plot <- ggplot() +
    # Plot aligned ordinations with semi-transparency
    geom_point(data = aligned_df, aes(x = Dim1, y = Dim2), color = "grey70", alpha = 0.3) +
    # Overlay the consensus ordination with distinct color and size
    geom_point(data = consensus_df, aes(x = Dim1, y = Dim2), color = "red", size = 2) +
    theme_minimal() +
    ggtitle("Aligned Ordinations with Consensus Overlaid") +
    xlab("Dimension 1") +
    ylab("Dimension 2")
  
  return(plot)
}

