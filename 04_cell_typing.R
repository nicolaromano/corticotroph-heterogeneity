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

plot_features_histo <- function(seurat_obj, feature, log_y_axis = TRUE) {
  expr <- data.frame(expr = GetAssayData(seurat_obj)[feature, ])

  if (log_y_axis) {
    ggplot(expr, aes(x = expr)) +
      geom_histogram(binwidth = 0.1) +
      scale_y_log10()
  } else {
    ggplot(expr, aes(x = expr)) +
      geom_histogram(binwidth = 0.1)
  }
}

type_cells <- function(seurat_obj) {
  # TODO add other cell types
}

pomc_hist <- lapply(seurat_objects, plot_features_histo, "Pomc")
marrangeGrob(pomc_hist, ncol = 4, nrow = 3)

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
