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

repeated_rarefaction <- function(physeq, repeats = 10, threshold = 250, method ="NMDS", colorb, shapeb, cloud = FALSE, ellipse = TRUE) {
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

  if (!(method == "NMDS")){
    stop("The only ordination method supported is NMDS.")
  }

  if (repeats <=4 & ellipse == TRUE){
    warning("Too few repeats to draw confidence ellipses.")
    ellipse <- F
  }

  sample_data(physeq)$sample_id <- rownames(sample_data(physeq))
  sample_id <- "sample_id"
  step1 <- rep_raref(as(otu_table(physeq), "matrix"), sample_data(physeq), threshold, repeats)
  step2 <- ord_and_mean(step1$repeat_count, step1$repeat_info, sample_data(physeq), repeats, method, sample_id)
  step3 <- plot_rep_raref(step2$ordinate_object, step2$physeq_object, step2$df_all, step2$df_median, colorb, shapeb, cloud, ellipse)
  print(step3)
  return(invisible(list("repeat_count" = step1$repeat_count, "repeat_info" = step1$repeat_info, "ordinate_object" = step2$ordinate_object, "physeq_object" = step2$physeq_object, "df_all" = step2$df_all, "df_median" = step2$df_median, "plot" = step3)))
}

#' Rarefaction is perforned reaptedly depending on input and a matrix containing all data is created together with a an info file reflecting it.
#'
#' @param count A matrix. An otu-table containing the count data the rarefaction should be performed on.
#' @param info Sample Data object from the phyloseq package.
#' Should contain info about the count data with matching column names.
#' @param threshold An integer. The threshold value at which to perform rarefaction.
#' @param repeats An integer. The amount of repeated rarefactions to perform.
#' A value = 1 means only one iteration of rarefaction is perfomed and therefore no repeats.
#' @returns A list containing a matrix with the repeated count table and another
#' matrix with the repeated info.
rep_raref <- function(count, info, threshold, repeats) {
  if (repeats == 0) {
    warning("repeats can't be 0. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  if (repeats < 0) {
    warning("repeats can't be negative. It needs to be a positive integer. Performs rarefaction without repetition.")
  }

  # Set up working files and perform one initial rarefaction
  rarified_count <- rrarefy(t(count), threshold)
  duplicated_info <- info
  sample_names_list1 <- row.names(info)
  sample_names_list2 <- colnames(count)

#  Perform repeated rarefaction
#  If repeats = 1, skips loop as no repetitions are needed.

  if (repeats > 1) {
    # Preallocate memory based on expected dimensions
    rarified_count <- matrix(NA, nrow = (repeats - 1) * nrow(rarified_count), ncol = ncol(rarified_count))
    duplicated_info <-matrix(NA, nrow = (repeats - 1) * nrow(info), ncol = ncol(info))
    sample_names_list1 <- character((repeats - 1) * nrow(info))
    sample_names_list2 <- character((repeats - 1) * ncol(count))

    # Initialize current row counters
    current_row_rarified <- 1
    current_row_info <- 1

    for (x in 2:repeats) {
      # Perform rarefaction
      rarified_result <- rrarefy(t(count), threshold)

      # Assign to preallocated matrix, ensuring range matches the number of rows
      rarified_count[current_row_rarified:(current_row_rarified + nrow(rarified_result) - 1),] <- rarified_result
      duplicated_info[current_row_info:(current_row_info + nrow(info) - 1),] <- as.matrix(info)

      # Append sample names for this iteration
      sample_names_list1[((x - 2) * nrow(info) + 1):((x - 1) * nrow(info))] <- paste0(row.names(info), "_", x)
      sample_names_list2[((x - 2) * ncol(count) + 1):((x - 1) * ncol(count))] <- paste0(colnames(count), "_", x)

      # Update the column indices
      current_row_rarified <- current_row_rarified + nrow(rarified_result)
      current_row_info <- current_row_info + nrow(info)
    }

    rownames(rarified_count) <- sample_names_list1
    colnames(rarified_count) <- rownames(count)
    rownames(duplicated_info) <- sample_names_list2
    colnames(duplicated_info) <- colnames(info)

    # Convert preallocated matrices to data frames if needed
    rarified_count <- as.data.frame(rarified_count)
    duplicated_info <- as.data.frame(duplicated_info)

  }

  rarified_count <- decostand(rarified_count, "normalize")

  # Transform both count data and info into dataframes
  # Corrects row names for compatibility in ordination function.
  duplicated_info <- as.data.frame(duplicated_info)
  row.names(duplicated_info) <- sample_names_list1
  rarified_count <- as.data.frame(rarified_count)
  row.names(rarified_count) <- sample_names_list2

  return(list("repeat_count" = rarified_count, "repeat_info" = duplicated_info))
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
ord_and_mean <- function(repeat_count, repeat_info, sample_info, repeats, method, sample_id) {
  if (!(sample_id %in% colnames(repeat_info))) {
    stop(paste0(sample_id, " is not a column name in input for repeat_info."))
  }

  # Convert the input into physeq objects
  rare_count_phy_repeat <- otu_table(t(repeat_count), taxa_are_rows = T)
  sample_info_tab_phy_repeat <- sample_data(repeat_info)
  rare_physeq_repeat <- phyloseq(rare_count_phy_repeat, sample_info_tab_phy_repeat)

  # Perform the Ordination calculation
  vst_pcoa <- phyloseq::ordinate(rare_physeq_repeat, method = method, distance = "bray")

  # Convert ordination result into a data frame object
  my_plot <- (plot_ordination(rare_physeq_repeat, vst_pcoa, justDF = TRUE))

  # Create matrix with median position for each sample
  group_by_sample <- my_plot %>%
    group_by(sample_id) %>%
    summarise(
      NMDS1 = median(NMDS1),
      NMDS2 = median(NMDS2)
    )

  df_median <- merge(group_by_sample, data.frame(sample_info), sample_id)

  return(list("ordinate_object" = vst_pcoa, "physeq_object" = rare_physeq_repeat, "df_all" = my_plot, "df_median" = df_median))
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
plot_rep_raref <- function(ordination, physeq, all_positions, median_positions, color, shape, cloud, ellipse) {
  color_median <- unlist(median_positions[color])
  color_tot <- unlist(all_positions[color])

  shape_median <- unlist(median_positions[shape])
  shape_tot <- unlist(all_positions[shape])
  
  # Create the plot and print it
  p <- plot_ordination(physeq, ordination, color = color, shape = shape)
  if (!cloud) { p$layers <- p$layers[-1] }
  p <- p + {if (ellipse) stat_ellipse(linetype = 1, lwd = 0.8, aes(color = color_tot, group = all_positions$sample_id))} +
    {if (!cloud) geom_point(data = median_positions, mapping = aes(x = median_positions[, 2], y = median_positions[, 3], color = color_median, shape = shape_median))}

  return(p)
}

