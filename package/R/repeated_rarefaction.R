repeated_rarefaction <- function(physeq, repeats, threshold, method, sample_id, colorb, shapeb, cloud = FALSE, ellipse = TRUE) {
  step1 <- rep_raref(as(otu_table(physeq), "matrix"), sample_data(physeq), threshold, repeats)
  step2 <- ord_and_mean(step1$repeat_count, step1$repeat_info, sample_data(physeq), repeats, method, sample_id)
  step3 <- plot_rep_raref(step2$ordinate_object, step2$physeq_object, step2$df_all, step2$df_mean, colorb, shapeb, cloud, ellipse)

  return(list("repeat_count" = step1$repeat_count, "repeat_info" = step1$repeat_info, "ordinate_object" = step2$ordinate_object, "physeq_object" = step2$physeq_object, "df_all" = step2$df_all, "df_median" = step2$df_median, "plot" = step3))
}

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

  # Perform repeated rarefaction
  # If repeats = 1, skips loop as no repetitions are needed.
  if (repeats > 1) {
    for (x in 2:repeats) {
      rarified_count <- rbind(rarified_count, rrarefy(t(count), threshold))
      duplicated_info <- rbind(duplicated_info, info)
      sample_names_list1 <- append(sample_names_list1, paste0(row.names(info), paste0("_", x)))
      sample_names_list2 <- append(sample_names_list2, paste0(colnames(count), paste0("_", x)))
    }
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

ord_and_mean <- function(repeat_count, repeat_info, sample_info, repeats, method, sample_id) {
  if (!(sample_id %in% colnames(repeat_info))) {
    stop(paste0(sample_id, " is not a column name in input for repeat_info."))
  }

  # Convert the input into physeq objects
  rare_count_phy_repeat <- otu_table(t(repeat_count), taxa_are_rows = T)
  sample_info_tab_phy_repeat <- sample_data(repeat_info)
  rare_physeq_repeat <- phyloseq(rare_count_phy_repeat, sample_info_tab_phy_repeat)

  # Perform the Ordination calculation
  vst_pcoa <- ordinate(rare_physeq_repeat, method = method, distance = "bray")

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
