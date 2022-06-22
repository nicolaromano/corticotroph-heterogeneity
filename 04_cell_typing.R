library(Seurat)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)

datasets <- read.csv("datasets.csv")

filenames <- dir("rds_outs/", pattern = "SCT.rds")
filenames <- paste0("rds_outs/", filenames)

# Read all datasets
# seurat_objects <- lapply(filenames, function(f)
#   {
#   print(paste("Reading", f))
#   seuratobj <- readRDS(f)
#   return(seuratobj)
#   })

# PCA elbow plots
plots <- lapply(seurat_objects, function(s) {
  p <- ElbowPlot(s, 30)
})
marrangeGrob(plots, nrow = 3, ncol = 4)

show_markers <- function(seurat_obj, hormonal = TRUE) {
  #' Plots [non]-hormonal features for the given dataset
  #' @param seurat_obj Seurat object
  #' @param hormonal Whether to plot features of hormonal cell
  #' types or those of other cell types
  #' @return A plot of the features

  if (hormonal) {
    features <- c(
      "Gh1", # Somatotrophs
      "Prl", # Lactotrophs
      # Corticotrophs (no Pax7 and no Pcsk2) + Melanotrophs
      "Pomc", "Pax7", "Pcsk2",
      "Lhb", "Fshb", # Gonadotrophs
      "Tshb" # Thyrotrophs
    )
  } else {
    features <- c(
      "Sox2", # Sox2 +ve cells
      "Pecam1", # Endothelial cells
      "Col1a1", # Fibroblasts
      "Cspg4", # Pericytes
      "Ptprc" # WBC
    )
  }

  FeaturePlot(seurat_obj, features, pt.size = 1) &
    scale_color_viridis_c() & # guide="none") &
    xlab(expression(UMAP[1])) &
    ylab(expression(UMAP[2]))
}

plot_features_histo <- function(seurat_obj, feature, find_threshold = TRUE,
                                log_y_axis = TRUE) {
  if (!feature %in% rownames(seurat_obj)) # Feature does not exist
    {
    p <- ggplot() +
      ggtitle(paste(seurat_obj$orig.ident[1], "- gene not found"))
    
    return(p)
    }
  
  expr <- data.frame(expr = GetAssayData(seurat_obj)[feature, ])

  p <- ggplot(expr, aes(x = expr)) +
    geom_histogram(binwidth = 0.1) +
    xlim(0, 10) +
    xlab(paste(feature, "expression")) +
    ggtitle(seurat_obj$orig.ident[1])
    
  if (log_y_axis) {
      p <- p + scale_y_log10()
  } 
  
  if (find_threshold)
    {
    thr <- otsu_thresh(expr$expr)
    p <- p + geom_vline(xintercept = thr, col = "red", lty = "dashed")
    }
  
  p
}

otsu_thresh <- function(values, levels = 100)
  {
  #' A simple implementation of Otsu's threshold
  #' Adapted from EBImage::otsu
  #' @param values: an array of values
  #' @param levels: the number of levels to consider
  #' @return : the Otsu's threshold for the values

  range <- c(0, max(values))
  breaks <- seq(range[1], range[2], length.out = levels + 1)
  
  h <- hist(values, breaks = breaks, plot = FALSE)
  counts <- as.double(h$counts)
  mids <- as.double(h$mids)
  len <- length(counts)
  w1 <- cumsum(counts)
  w2 <- w1[len] + counts - w1
  cm <- counts * mids
  m1 <- cumsum(cm)
  m2 <- m1[len] + cm - m1
  var <- w1 * w2 * (m2/w2 - m1/w1)^2
  maxi <- which(var == max(var, na.rm = TRUE))
  
  return((mids[maxi[1]] + mids[maxi[length(maxi)]])/2)
  }
  
type_cells_thr <- function(seurat_obj) {
  # TODO add other cell types
  # We call corticotrophs cells that 
}

pomc_hist <- lapply(seurat_objects, plot_features_histo, "Pomc", 
                    find_threshold = TRUE, log_y_axis = TRUE)
do.call("grid.arrange", c(pomc_hist, ncol = 4, nrow = 3))

pax7_hist <- lapply(seurat_objects, plot_features_histo, "Pax7", 
                    find_threshold = TRUE, log_y_axis = TRUE)
do.call("grid.arrange", c(pax7_hist, ncol = 4, nrow = 3))

seurat_obj <- seurat_objects[[3]]
# Try different resolutions and use the elbow method to find the "optimal one"
seurat_obj %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.9) -> seurat_obj
DimPlot(seurat_obj, label = TRUE, pt.size = 1)

show_markers(seurat_obj, hormonal = TRUE)
show_markers(seurat_obj, hormonal = FALSE)

FeaturePlot(seurat_obj, features_other, pt.size = .1) &
  scale_color_viridis_c() &
  xlab(expression(UMAP[1])) &
  ylab(expression(UMAP[2]))
