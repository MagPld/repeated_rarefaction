library(vegan)
library(geomorph)
#library(ggpubr)
library(doParallel)
library(foreach)
library(ggplot2)

rarefaction_depth <- 400 
# Initialize a list to store the rarefied matrices
rarefied_matrices <- list()
count_matrix <- data.frame(t(.otu_table))

# Perform rarefaction n times and store each result
for (i in 1:400) {
  rarefied_matrices[[i]] <- rrarefy(count_matrix, sample = rarefaction_depth)
}

#========================= ordinations and plots generation

# Initialize a list to store ordination plots
ordination_plots <- list()
ordinations <- list()

# Set up parallel backend to use available cores
numCores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Perform parallel computation using foreach
results <- foreach(i = 1:length(rarefied_matrices), .packages = c('vegan', 'ggplot2')) %dopar% {
  # Calculate Bray-Curtis distance
  dist_matrix <- vegdist(rarefied_matrices[[i]], method = "bray")
  
  # Perform PCoA (Principal Coordinates Analysis)
  ordination <- cmdscale(dist_matrix, k = 2)
  
  # Convert ordination result to data frame
  ord_df <- as.data.frame(ordination)
  colnames(ord_df) <- c("PCoA1", "PCoA2")
  ord_df$Sample <- rownames(ordination)
  
  # Create a ggplot for the ordination
  p <- ggplot(ord_df, aes(x = PCoA1, y = PCoA2)) +
    geom_point() +
    ggtitle(paste("Ordination", i)) +
    theme_minimal()
  
  # Return a list containing the ordination and the plot
  list(ordination = ordination, plot = p)
  
}

# Stop the cluster after computation
stopCluster(cl)


# Extract ordinations and plots from the results
for (i in 1:length(results)) {
  ordinations[[i]] <- results[[i]]$ordination
  ordination_plots[[i]] <- results[[i]]$plot
}

#================================ procrustes

# Perform Procrustes analysis to align all ordinations to the first one
aligned_ordinations <- lapply(ordinations, function(x) procrustes(ordinations[[1]], x)$Yrot)

# Convert list to array for consensus calculation
aligned_array <- array(unlist(aligned_ordinations), dim = c(nrow(aligned_ordinations[[1]]), ncol(aligned_ordinations[[1]]), length(aligned_ordinations)))

# Compute consensus using mean shape
consensus_coords <- mshape(aligned_array)

# Plot the consensus ordination
plot(consensus_coords)

#=============================== plotting

# Combine aligned ordinations into one data frame for plotting
aligned_df <- data.frame()

for (i in 1:length(aligned_ordinations)) {
  temp_df <- as.data.frame(aligned_ordinations[[i]])
  colnames(temp_df) <- c("Dim1", "Dim2")
  temp_df$Sample <- rownames(aligned_ordinations[[i]])
  temp_df$Ordination <- paste0("Ordination", i)
  aligned_df <- rbind(aligned_df, temp_df)
}

# Convert consensus_coords to a data frame
consensus_df <- as.data.frame(consensus_coords)
colnames(consensus_df) <- c("Dim1", "Dim2")
consensus_df$Sample <- rownames(consensus_coords)

# Plot all aligned ordinations
ggplot() +
  # Plot aligned ordinations with semi-transparency
  geom_point(data = aligned_df, aes(x = Dim1, y = Dim2), color = "grey70", alpha = 0.3) +
  # Overlay the consensus ordination with distinct color and size
  geom_point(data = consensus_df, aes(x = Dim1, y = Dim2), color = "red", size = 2) +
  theme_minimal() +
  ggtitle("Aligned Ordinations with Consensus Overlaid") +
  xlab("Dimension 1") +
  ylab("Dimension 2")
