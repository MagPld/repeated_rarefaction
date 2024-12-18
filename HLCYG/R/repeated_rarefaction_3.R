#' Perform one iteration of repeated rarefaction and produce ordination plot
#'
#' @param physeq A phyloseq object.
#' @param repeats A positive integer. Indicates the amount of repeats run.
#' A value = 1 means no repeats are used. Default = 50.
#' @param threshold A positive integer. The threshold value for the rarefaction. Default = 250
#' @param colorb A string. Name of the column with data to colour the graph by.
#' @param shapeb A string. Name of the column with data to shape the points on the graph by.
#' @param cloud Boolean. Aesthetic setting for the graph. Default is FALSE.
#' TRUE shows datapoints for all repeats.
#' @param ellipse Boolean. Aesthetic setting for the graph. Default is TRUE.
#' @param cores An int. The amount of cores used when running. The default = 4.
#' @returns A list with repeat count table, repeat info table, ordination object,
#' physeq object, dataframe with all repeat ordination positions, dataframe with median ordination positions, and ordination plot.
#' @examples
#' repeated_rarefaction(HLCYG_physeq_data, repeats=5, threshold=500, method="NMDS", colorb="sample_id", shapeb="location", T, F)
#' repeated_rarefaction(HLCYG_physeq_data, repeats=10, threshold=250, method="NMDS", colorb="sample_id", shapeb="location", T, T)


repeated_rarefaction_2 <- function(input, repeats = 50, threshold = 250, colorb="sample_id", group="sample_id", cloud = TRUE, ellipse = FALSE, cores = 4) {
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

  # ============ Checks and warnings

  if (!(colorb %in% names(sample_data(physeq)))) {
    stop(paste("'",colorb,"' is not a column name in the sample information in the inputed phyloseq object.
                  repeated_rarefaction needs an existing column to color the ordination plot by.", sep=""))
  }
  if (!(group %in% names(sample_data(physeq)))) {
    stop(paste("'",group,"' is not a column name in the sample information in the inputed phyloseq object.
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
  step1 <- rep_raref_2(data.frame(t(otu_table(physeq))), threshold, repeats)
  step2 <- ord_and_mean_2(step1$rarefied_matrix_list, repeats, cores)
  step3 <- plot_rep_raref_2(step2$aligned_ordinations, step2$consensus_coordinates, sample_data(physeq), colorb, group, cloud, ellipse, "Aligned Ordinations with Consensus Overlaid")

  print(step3$plot)

  return(invisible(list("repeats" = repeats, "df_consensus_coordinates" = step3$consensus_df, "df_all" = step3$df_all, "plot" = step3$plot)))
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
rep_raref_2 <- function(count, threshold, repeats) {
  if (repeats == 0) {
    warning("repeats can't be 0. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  if (repeats < 0) {
    warning("repeats can't be negative. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  # Set up working files
  rarefied_matrices <- list()

  # Perform repeated rarefaction and store the normalized results in a list
  if (repeats >= 1) {
    for (i in 1:repeats) {
      rarefied_count <- rrarefy(count, sample = threshold)
      rarefied_matrices[[i]] <- rarefied_count
    }
  }

  return(invisible(list("rarefied_matrix_list"=rarefied_matrices)))
}

#' Calculates ordination on the repeated count input and as well as a median result for each original datapoint.
#'
#' @param rarified_matrix_list A matrix. Contains the repeated count data.
#' @param repeats An integer. The amount of repeats that has been performed on the data.
#' @returns Returns a list containg an ordination object, a phyloseq object with
#' the repeated count data, a dataframe with all positions from the ordination
#' calculaction, and a datafram with just the median position from the calculation.
ord_and_mean_2 <- function(rarefied_matrix_list, repeats, cores = 4) {

  #========================= ordinations and plots generation

  # Initialize a list to store ordinations
  ordinations <- list()

  # Set up parallel backend to use available cores
  # numCores <- detectCores() - 1  , Use one less than the total number of cores
  cl <- makeCluster(cores)
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

  return(invisible(list("aligned_ordinations" = aligned_ordinations, "consensus_coordinates" = consensus_coords)))
}


#' Create a plot for the repeated rarefaction.
#'
#' @param ordination An ordination object.
#' @param physeq A phyloseq object containing the data to be plotted.
#' @param all_positions A dataframe containing position data for all repeats.
#' @param mediant_positions A dataframe containing median position data for each original data point.
#' @param color A string. Name of the column to color points by.
#' @param shape A string. Name of the column to shape points by. Substitute by "group"
#' @param cloud Boolean. Aesthetic setting for the graph. Default is FALSE.
#' TRUE shows datapoints for all repeats.
#' #' @param ellipse Boolean. Aesthetic setting for the graph. Default is TRUE.
#' @returns Returns an ordination plot.
plot_rep_raref_2 <- function(aligned_ordinations, consensus_coordinates, info, color, group, cloud, ellipse, title) {

  # Combine aligned ordinations into one data frame for plotting
  aligned_df <- data.frame()

  for (i in 1:length(aligned_ordinations)) {
    temp_df <- as.data.frame(aligned_ordinations[[i]])
    colnames(temp_df) <- c("Dim1", "Dim2")
    temp_df$sample_id <- rownames(aligned_ordinations[[i]])
    temp_df$ordination <- paste0("Ordination", i)
    temp_df[[color]] <- info[[color]]
    temp_df[[group]] <- info[[group]]
    aligned_df <- rbind(aligned_df, temp_df)
  }

  # Convert consensus_coordinates to a data frame
  consensus_df <- as.data.frame(consensus_coordinates)
  colnames(consensus_df) <- c("Dim1", "Dim2")
  consensus_df$sample_id <- rownames(aligned_ordinations[[1]])

  # Determine how many colors are needed to plot the levels of the color variable
  num_levels <- length(unique(consensus_df[[color]]))

  plot <- ggplot()

  if (num_levels <= 6) {
    # Use variable-based coloring if there are 6 or fewer categories
    plot <- plot +
      geom_point(
        data = aligned_df,
        aes(x = Dim1, y = Dim2, color = .data[[color]]),
        alpha = 0.3
      )
  } else {
    # Use a fixed grey color for all points if there are more than 6 categories
    plot <- plot +
      geom_point(
        data = aligned_df,
        aes(x = Dim1, y = Dim2),
        color = "grey70",
        alpha = 0.3
      )
    }

  if (!cloud) {
    plot$layers <- plot$layers[-1]
  }


  if (ellipse) {
    if (num_levels <= 6) {
      # Ellipses colored by the color variable
      plot <- plot +
        stat_ellipse(
          data = aligned_df,
          aes(x = Dim1, y = Dim2, color = .data[[color]], group = .data[[group]]),
          linetype = 1, lwd = 0.8
        )
    } else {
      # Ellipses in grey if more than 6 categories
      # Still grouping by the color variable to get multiple ellipses if there are multiple groups
      plot <- plot +
        stat_ellipse(
          data = aligned_df,
          aes(x = Dim1, y = Dim2, group = .data[[group]]),
          color = "grey70",
          linetype = 1, lwd = 0.8
        )
    }
  }

  # Plot all aligned ordinations (consensus points)
  plot <- plot +
    geom_point(data = consensus_df, aes(x = Dim1, y = Dim2), color = "red", size = 2) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(title) +
    xlab("Dimension 1") +
    ylab("Dimension 2")

  return(invisible(list("plot" = plot, "consensus_df" = consensus_df, "df_all" = aligned_df)))
}



