.onLoad <- function(libname, pkgname) {
    bioc_packages <- c("phyloseq") # List of Bioconductor packages
    missing_packages <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
    
    if (length(missing_packages) > 0) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(missing_packages)
    }
    
    # Disable "rgl" creation of graphical windows, since we do not require
    # this functionality from "geomorph".
    options(rgl.useNULL = TRUE)
    
    library(phyloseq)
    library(clusterSim)
    library(dplyr)
    library(ggplot2)
    library(vegan)
    library(geomorph)
    library(doParallel)
    library(foreach)
}
