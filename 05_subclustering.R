library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pbapply)
library(pheatmap)
library(parallel)

#### SUBSET DATA TO COMMON GENES ####

# First we load all data (M and F) and find the common
# genes between all datasets. We then subset all datasets to these genes.
# This will save some headaches later on.
filenames <- dir("rds_outs", pattern = "corticotrophs.rds", full.names = TRUE)
# Remove Ho2020, which has only a handful of cells.
# This also increases the number of common genes.
filenames <- filenames[-grep("Ho2020", filenames)]
all_corticotrophs <- pblapply(filenames, readRDS)
names(all_corticotrophs) <- filenames

gene_names <- lapply(all_corticotrophs, function(x) {
  rownames(x)
})
# 11103 common genes
common_genes <- Reduce(intersect, gene_names)

all_corticotrophs <- lapply(all_corticotrophs, function(x) {
  x <- x[common_genes, ]
})

# Save the subset files
for (i in 1:length(all_corticotrophs)) {
  print(paste("Saving", names(all_corticotrophs)[i]))
  saveRDS(all_corticotrophs[[i]], names(all_corticotrophs)[i])
}

# Sanity-check: do all datasets have the same number of genes?
sapply(all_corticotrophs, function(x) {
  paste(nrow(x), "genes x", ncol(x), "cells")
})

# We can remove the data from memory now
rm(all_corticotrophs)
gc()

#### SUBCLUSTERING ####

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Data to process ("M" or "F")
data_to_process <- "M"

plot_out_type <- "none" # "png", "pdf" or "none"

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

names(seurat_corticotrophs) <- sapply(seurat_corticotrophs, function(x) {
  x$orig.ident[1]
})

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
  print(paste("Elbow method for", obj$author[1], obj$year[1], "-", obj$sex[1]))
  print("***************************************************")

  wss <- lapply(resolutions, function(resol) {
    print("***************************************")
    print(paste("Clustering at resolution", resol))
    # The results of this clustering get stored in:
    # 1) obj$seurat_clusters (which is a shortcut for
    #    obj@metadata$seurat_clusters)
    # 2) obj$SCT_snn_res.xxx where xxx is the resolution
    # Option 1) is more convenient to use in this function,
    # but mind that it gets overwritten every time FindClusters is run!
    # Note we use mclapply to parallelize the process
    # We use all available cores (detectCores())
    res <- mclapply(1:n_shuffle,
      mc.cores = detectCores(),
      FUN = function(x) {
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
      }
    )

    # Append to results
    wss <- data.frame(
      dataset = paste(obj$author[1], obj$year[1], "-", obj$sex[1]),
      resol = resol,
      wss = sapply(res, function(x) {
        mean(x$wss)
      }),
      min_n_clusters = min(sapply(res, nrow)),
      max_n_clusters = max(sapply(res, nrow))
    )

    wss
  })

  do.call(rbind, wss) %>%
    as.data.frame() %>%
    pivot_longer(starts_with("wss"), names_to = "shuffle", values_to = "wss") %>%
    mutate(shuffle = as.numeric(str_extract(shuffle, "[0-9]+"))) %>%
    group_by(dataset, resol) %>%
    summarise(
      wss_mean = mean(wss),
      wss_sd = sd(wss),
      min_n_clusters = min(min_n_clusters),
      max_n_clusters = max(max_n_clusters)
    ) %>%
    ungroup()
}

resolutions <- seq(0.1, 1.5, 0.1)

if (calculate_elbow_plots) {
  # This is quite slow (~7 minutes using 16 cores)
  start_time <- Sys.time()
  elbow_plot_data <- lapply(seurat_corticotrophs, function(obj) {
    obj <- obj %>%
      RunUMAP(dims = 1:20) %>%
      elbow_method(
        resolutions = resolutions,
        n_shuffle = 100, subsampling = 0.6
      )
  })
  end_time <- Sys.time()
  print(end_time - start_time)

  # Save results to file
  write.csv(do.call(rbind, elbow_plot_data),
    file = paste0("elbow_plots_", data_to_process, ".csv"),
    row.names = FALSE
  )
}

wss_data <- read.csv(paste0(
  "elbow_plots_",
  data_to_process, ".csv"
))

if (plot_out_type == "png") {
  png(paste0("plots/elbow_plots_", data_to_process, ".png"),
    width = 15, height = 10, units = "in", res = 300
  )
} else if (plot_out_type == "pdf") {
  pdf(paste0("plots/elbow_plots_", data_to_process, ".pdf"),
    width = 15, height = 10
  )
}

# Plot the results
ggplot(wss_data, aes(x = resol, y = wss_mean)) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = wss_mean - wss_sd,
      ymax = wss_mean + wss_sd
    ),
    alpha = 0.2
  ) +
  geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max(resolutions))) +
  geom_text(aes(label = paste(min_n_clusters, max_n_clusters, sep = "-")), vjust = -1.5,
  size = 4) +
  facet_wrap(~dataset, ncol = 3, scales = "free") +
  xlab("Resolution") +
  ylab("WSS") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

if (plot_out_type != "none") {
  dev.off()
}

# These were manually chosen based on the elbow plots
# These are the values where the WSS change is < 10%
if (data_to_process == "M") {
  resolutions <- c(0.5, 0.4, 0.6, 0.5, 0.5, 0.4)
} else if (data_to_process == "F") {
  resolutions <- c(0.4, 0.3, 0.4, 0.5)
}

seurat_corticotrophs <- lapply(
  seq_along(seurat_corticotrophs),
  function(i) {
    print(paste0(
      "Clustering ", seurat_corticotrophs[[i]]$orig.ident[1],
      " at resolution ", resolutions[i]
    ))
    obj <- seurat_corticotrophs[[i]]
    obj <- obj %>%
      RunUMAP(dims = 1:20) %>%
      FindNeighbors(k = 30) %>%
      FindClusters(resolution = resolutions[i])

    return(obj)
  }
)

cort_plots <- lapply(seurat_corticotrophs, function(obj) {
  DimPlot(obj, pt.size = 1.5) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(paste(obj$author[1], obj$year[1], "-", obj$sex[1])) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      title = element_text(size = 16)
    )
})

if (plot_out_type == "png") {
  png(paste0("plots/corticotrophs_all_datasets_", data_to_process, ".png"),
    width = 15, height = 10, units = "in", res = 300
  )
} else if (plot_out_type == "pdf") {
  pdf(paste0("plots/corticotrophs_all_datasets_", data_to_process, ".pdf"),
    width = 15, height = 10
  )
}
do.call("grid.arrange", c(cort_plots, ncol = 3))

if (plot_out_type != "none") {
  dev.off()
}

# Save to RDS
pblapply(seurat_corticotrophs, function(obj) {
  saveRDS(obj, file = paste0("rds_outs/", obj$orig.ident[1], "_subclustered.rds"))
})
