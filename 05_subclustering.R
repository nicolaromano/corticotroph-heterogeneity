library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pbapply)
library(pheatmap)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Data to process ("M" or "F")
data_to_process <- "F"

# This is very heavy to compute, so this is a switch to turn it off if already computed
calculate_elbow_plots <- FALSE

# Read only M or F files
filenames <- dir("rds_outs",
  pattern = paste0(data_to_process, "_corticotrophs.rds"),
  full.names = TRUE
)

# Ho 2020 has only a handful of cells - we'are not considering it
filenames <- filenames[-grep("Ho2020", filenames)]
seurat_corticotrophs <- pblapply(filenames, readRDS)

prettify_df_name <- function(name, show_sex = TRUE) {
  name <- as.character(name)
  paste(
    substr(name, 1, nchar(name) - 5), # Author
    substr(name, nchar(name) - 4, nchar(name)-1),
    ifelse (show_sex,
      substr(name, nchar(name), nchar(name)),
      "")
  ) # Year
}

##### CLUSTERING QUALITY ASSESSMENT #####

### Elbow method ###
# We look at the "elbow" in a plot of resolution vs WSS
# (Within-Cluster-Sum of Squared Errors)
elbow_method <- function(obj, resolutions = seq(0, 1, 0.1),
                         n_shuffle = 100, subsampling = 0.8) {
  #' Calculated the WSS at various clustering resolutions
  #' This performs clustering on the data at various resolutions
  #' and calculates the WSS (Within-Cluster-Sum of Squared Errors)
  #' for each resolution. The process is performed on random subset of the data.
  #' @param obj: the Seurat object
  #' @param resolutions: the resolutions to test (default is seq(0, 1, 0.1))
  #' @param n_shuffle: the number of times to shuffle the data
  #' @param subsampling: the fraction of cells to subsample each time
  #' @return a data frame with the resolution and WSS for each resolution
  
  print("***************************************************")
  print(paste0("Elbow method for ", prettify_df_name(obj$orig.ident[1])))
  print("***************************************************")

  wss <- lapply(resolutions, function(resol){
    print("***************************************")
    print(paste("Clustering at resolution", resol))
    # The results of this clustering get stored in:
    # 1) obj$seurat_clusters (which is a shortcut for
    #    obj@metadata$seurat_clusters)
    # 2) obj$SCT_snn_res.xxx where xxx is the resolution
    # Option 1) is more convenient to use in this function,
    # but mind that it gets overwritten every time FindClusters is run!
    res <- lapply(1:n_shuffle, function(x) {
      # Subsample the data - Note the use of seed = NULL
      # to ensure that the subsampling is random every time
      # See https://github.com/mojaveazure/seurat-object/issues/62 for details
      print(paste("ROUND", x))
      sub_obj <- obj[, sample(Cells(obj),
        size = subsampling * ncol(obj),
        replace = FALSE
      ), seed = NULL]
      # Run clustering
      sub_obj <- sub_obj %>% 
        FindNeighbors() %>% 
        FindClusters(resolution = resol)

      # Find the WSS
      res <- data.frame(center = colMeans(GetAssayData(sub_obj))) %>%
        rownames_to_column("cell") %>%
        # Add the cluster ID
        mutate(cluster = sub_obj$seurat_clusters) %>%
        # Calculate the WSS for each cluster
        group_by(cluster) %>%
        summarise(wss = sum((center - mean(center))^2)) %>%
        ungroup()
    })

    # Append to results
    wss <- data.frame(
      resol = resol,
      wss = sapply(res, function(x) {
        mean(x$wss)
      }),
      min_n_clusters = min(sapply(res, nrow)),
      max_n_clusters = max(sapply(res, nrow))
    )

    wss
  })

  do.call (rbind, wss) %>% 
    as.data.frame()  %>% 
    pivot_longer(starts_with("wss"), names_to = "shuffle", values_to = "wss")  %>% 
    mutate(shuffle = as.numeric(str_extract(shuffle, "[0-9]+"))) %>% 
    group_by(resol) %>% 
    summarise(wss_mean = mean(wss), 
      wss_sd = sd(wss),
      min_n_clusters = min(min_n_clusters), 
      max_n_clusters = max(max_n_clusters)) %>%
    ungroup() %>% 
    ggplot(aes(x = resol, y = wss_mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = wss_mean - wss_sd, ymax = wss_mean + wss_sd),
      alpha = 0.2) + 
    geom_point(size = 2) +
    geom_text(aes(label = paste(min_n_clusters, max_n_clusters, sep="-")), vjust = -1) +
    xlab("Resolution") +
    ylab("WSS") +
    # Works best with uniformly spaced (and not too many) resolutions
    scale_x_continuous(breaks = seq(min(resolutions), max(resolutions),
      length.out = length(resolutions)
    )) +
    ggtitle(prettify_df_name(obj$orig.ident[1], show_sex = TRUE)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 16)
    )
}

if(calculate_elbow_plots)
  {
  # This is very slow (~50 minutes) 
  # Could be parallelised but mclapply fails for some obscure reason...
  start_time <- Sys.time()
  elbow_plots <- lapply(seurat_corticotrophs, function(obj) {
    obj <- obj %>%
      RunUMAP(dims = 1:20) %>%
      elbow_method(resolutions = seq(0.1, 1.5, 0.1),
      n_shuffle = 100, subsampling = 0.6)
  })
  end_time <- Sys.time()
  print(end_time - start_time)

  do.call("grid.arrange", c(elbow_plots, ncol = 3))
    

  png(paste0("plots/elbow_plots", data_to_process, ".png"), width = 1600, height = 1200)
  do.call("grid.arrange", c(elbow_plots, ncol = 3))
  dev.off()
  }

# These were manually chosen based on the elbow plots
if data_to_process == "M"
  resolutions <- c(0.4, 0.3, 0.5, 0.3, 0.5, 0.4)
else if data_to_process == "F"
  resolutions <- c(0.4, 0.3, 0.4, 0.5)

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

do.call("grid.arrange", cort_plots)

png(paste0("plots/corticotrophs_all_datasets_", data_to_process, ".png"), 
  width = 800, height = 1200)
do.call("grid.arrange", cort_plots)
dev.off()

# Save to RDS
pblapply(seurat_corticotrophs, function(obj) {
  saveRDS(obj, file = paste0("rds_outs/", obj$orig.ident[1], "_subclustered.rds"))
})