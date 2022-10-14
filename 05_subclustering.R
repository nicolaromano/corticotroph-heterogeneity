library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pbapply)
library(pheatmap)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Read only male files
filenames <- dir("rds_outs",
  pattern = "M_corticotrophs.rds",
  full.names = TRUE
)

# Ho 2020 has only a handful of cells - we'are not considering it
filenames <- filenames[-grep("Ho2020", filenames)]
seurat_corticotrophs <- pblapply(filenames, readRDS)

##### CLUSTERING QUALITY ASSESSMENT #####

### Elbow method ###
# We look at the "elbow" in a plot of resolution vs WSS
# (Within-Cluster-Sum of Squared Errors)
elbow_method <- function(obj, 
  resolutions = seq(0, 1, 0.1) {
  #' Calculated the WSS at various clustering resolutions
  #' This performs clustering on the data at various resolutions
  #' @param obj: the Seurat object
  #' @param resolutions: the resolutions to test (default is seq(0, 2, 0.1))
  wss <- NULL
  for (resol in resolutions)
  {
    print(paste("Clustering at resolution", resol))
    # The results of this clustering get stored in:
    # 1) obj$seurat_clusters (which is a shortcut for 
    #    obj@metadata$seurat_clusters)
    # 2) obj$SCT_snn_res.xxx where xxx is the resolution
    # Option 1) is more convenient to use in this function, 
    # but mind that it gets overwritten every time FindClusters is run!
    obj <- FindClusters(obj, resolution = resol)

    # Find the WSS
    res <- data.frame(center = colMeans(GetAssayData(obj))) %>% 
      rownames_to_column("cell") %>% 
      # Add the cluster ID
      mutate(cluster = obj$seurat_clusters) %>%
      # Calculate the WSS for each cluster
      group_by(cluster)  %>% 
      summarise(wss = sum((center - mean(center))^2)) %>%
      ungroup()

    # Append to results
    wss <- rbind(wss, c(
      resol = resol,
      wss = mean(res$wss),
      n_clusters = length(levels(obj$seurat_clusters))
    ))
  }

  wss <- data.frame(wss)

  ggplot(wss, aes(x = resol, y = wss)) +
    geom_line() +
    geom_point() +
    geom_text(aes(label = n_clusters), vjust = -1) +
    xlab("Resolution") +
    # Works best with uniformly spaced (and not too many) resolutions
    scale_x_continuous(breaks = seq(min(resolutions), max(resolutions),
      length.out = length(resolutions)
    )) +
    ggtitle(obj$orig.ident[1]) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 90)
    )
}

elbow_plots <- lapply(seurat_corticotrophs, function(obj) {
  obj <- obj %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors()

  elbow_method(obj, resolutions = seq(0.1, 1.5, 0.1))
})

do.call("grid.arrange", c(elbow_plots, ncol = 3))

# These were manually chosen
resolutions <- c(0.3, 0.3, 0.5, 0.5, 0.3, 0.3)

seurat_corticotrophs <- lapply(seq_along(seurat_corticotrophs), function(i) {
  obj <- seurat_corticotrophs[[i]]
  obj <- obj %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors() %>%
    FindClusters(resolution = resolutions[i])

  return(obj)
})

cort_plots <- lapply(seurat_corticotrophs, function(obj) {
  DimPlot(obj) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(obj$orig.ident[1])
})

png("plots/corticotrophs_all_datasets.png", width = 800, height = 1200)
do.call("grid.arrange", cort_plots)
dev.off()

# Save to RDS
pblapply(seurat_corticotrophs, function(obj) {
  saveRDS(obj, file = paste0("rds_outs/", obj$orig.ident[1], "_subclustered.rds"))
})
