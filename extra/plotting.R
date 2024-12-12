chunk_size <- 10
n_plots <- length(result$ordination_plots)
plot_grids <- list()
for (i in seq(1, n_plots, by = chunk_size)) {
  # Determine the end index of the current chunk
  end_idx <- min(i + chunk_size - 1, n_plots)
  
  # Extract the subset of plots for the current chunk
  current_plots <- result$ordination_plots[i:end_idx]
  
  # Wrap and display the grid of current chunk
  plot_grid <- wrap_plots(current_plots, ncol = 3)
  plot_grids[[length(plot_grids) + 1]] <- plot_grid
}

for (i in seq_along(plot_grids)) {
  filename <- paste0("plot_grid_", i, ".png")
  filepath <- paste0("/HLCYG/", filename)
  ggsave(filepath, plot = plot_grids[[i]], width = 10, height = 8, dpi = 300)
}
