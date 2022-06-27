library(Seurat)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(pbapply)

#### Various functions ####
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
      ggtitle(seurat_obj$orig.ident[1], subtitle = "Not found")
    
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
    thr <- otsu_thresh(expr$expr, levels = 500)
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
  # We call corticotrophs cells that are above the threshold for POMC
  expr <- GetAssayData(seurat_obj)["Pomc", ]
  thr <- otsu_thresh(expr)

  cells_bc <- Cells(subset(seurat_obj, Pomc > thr))  
  print(paste0(round(length(cells_bc)/length(Cells(seurat_obj)) * 100, 2), 
              "% cells typed as corticotrophs/melanotrophs"))
  levels(seurat_obj@active.ident) <- c(levels(seurat_obj@active.ident), "Cortico/Melano")
  Idents(seurat_obj)[cells_bc] <- "Cortico/Melano"
  p <- DimPlot(seurat_obj, pt.size = 0.5) +
    scale_color_manual(values = c("lightgray", "navy")) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    theme(legend.position = "none") +
    ggtitle(seurat_obj$orig.ident[1])
  
  return(p)
}

get_pomc_cells <- function(seurat_obj, resolution = 0.6, 
                            pomc_clusters = NULL, save_rds = TRUE)
  {
  seurat_obj %>%
    FindNeighbors() %>%
    FindClusters(resolution = resolution) -> seurat_obj
  
  p1 <- DimPlot(seurat_obj, label = TRUE, pt.size = .1) + 
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(seurat_obj$orig.ident[1]) +
    theme(legend.position = "none")
    
  p2 <- FeaturePlot(seurat_obj, feature = "Pomc", pt.size = .1) +
    scale_color_viridis_c() +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) 
  
  grid.arrange(p1, p2, nrow=2)
  
  if (save_rds)
    {
    seurat_obj %>% 
      subset(idents = pomc_clusters) %>% 
      RunUMAP(dims = 1:20) -> seurat_obj_pomc
    
    saveRDS(seurat_obj_pomc, file = paste0("rds_outs/",
                              seurat_obj_pomc$orig.ident[1],
                              "_pomc.rds"))
    return(seurat_obj_pomc)
    }
  
  return(seurat_obj)
  }

datasets <- read.csv("datasets.csv")

filenames <- dir("rds_outs/", pattern = "SCT.rds")
filenames <- paste0("rds_outs/", filenames)

# Read all datasets
seurat_objects <- pblapply(filenames, function(f)
{
  seuratobj <- readRDS(f)
  return(seuratobj)
})

# PCA elbow plots
plots <- lapply(seurat_objects, function(s) {
  p <- ElbowPlot(s, 30)
})

do.call("grid.arrange", c(plots, nrow = 3, ncol = 4))

#### POMC HISTOGRAM #####
pomc_hist <- lapply(seurat_objects, plot_features_histo, "Pomc", 
                    find_threshold = TRUE, log_y_axis = TRUE)
do.call("grid.arrange", c(pomc_hist, ncol = 4, nrow = 3))


pomc_cells <- lapply(seurat_objects, type_cells_thr)
do.call("grid.arrange", c(pomc_cells, ncol = 4, nrow = 3))

#### CLUSTERING ####

# These are manually picked as not possible to use the same for everyone
resolutions <- c(0.5, # Allensworth 2021
                 0.5, # Cheung 2018
                 0.3, 0.3, # Fletcher 2019 F/M
                 1, 1, # Ho 2020 F/M
                 0.5, # Kucka 2021
                 0.5, # Lopez 2021
                 0.5, # Mayran 2019
                 0.3, 0.3, # Ruf Zamoiski 2019 F/M
                 0.5) # Vennekens 2021

pomc_clusters <- list(8:9, # Allensworth 2021
                      c(6, 13), # Cheung 2018
                      3, c(3,9), # Fletcher 2019 F/M
                      c(4, 6, 8), c(4, 6, 7), # Ho 2020 F/M
                      c(7, 22, 23), # Kucka 2021
                      c(6, 9), # Lopez 2021
                      c(4, 5, 12), # Mayran 2019
                      c(4, 5, 19), c(2, 4, 17), # Ruf Zamoiski 2019 F/M
                      5) # Vennekens 2021

seurat_objects_pomc <- sapply(seq_along(seurat_objects), function(i) {
  seurat_obj <- seurat_objects[[i]]

  get_pomc_cells(seurat_obj, resolution = resolutions[i], 
                         pomc_clusters = pomc_clusters[[i]],
                 save_rds = TRUE)
  })

#### SPLIT POMC CELLS INTO CORTICOTROPHS AND MELANOTROPHS ####

# Load files (to start directly from here)
filenames <- dir(path = "rds_outs/", pattern = "pomc.rds")
seurat_objects_pomc <- pbsapply(paste0("rds_outs/", filenames), readRDS)

pomc_plots <- lapply(seurat_objects_pomc, function(seurat_obj) {
  DimPlot(seurat_obj, label = TRUE) +
    ggtitle(seurat_obj$orig.ident[1]) +
    theme(legend.position = "none")
})

do.call("grid.arrange", c(pomc_plots, ncol=4))

melano_plots <- lapply(seurat_objects_pomc, function(seurat_obj) {
  melano_markers <- c("Pcsk2", "Pax7")
  melano_markers <- melano_markers[melano_markers %in% rownames(seurat_obj)]
  
  if (length(melano_markers) > 1)
    melano_mark_expr <- data.frame(t(as.matrix(GetAssayData(seurat_obj)[melano_markers,])))
  else
    melano_mark_expr <- data.frame(as.matrix(GetAssayData(seurat_obj)[melano_markers,]))
  umap_coord <- data.frame(seurat_obj@reductions$umap@cell.embeddings)
  p <- ggplot(umap_coord, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(col = rowSums(melano_mark_expr) > 0), size = 0.5) +
    scale_color_manual(values = c("lightgray", "navy"), guide = "none") +
    ggtitle(seurat_obj$orig.ident[1])
  return(p)
  })

do.call("grid.arrange", c(melano_plots, nrow=3))

threshold <- 0.1

melanotrophs <- pbsapply(seurat_objects_pomc, function(seurat_obj)
  {
  # Pax7 is not always detected
  if ("Pax7" %in% rownames(seurat_obj))
    melano <- subset(seurat_obj, Pcsk2 >= threshold | Pax7 >= threshold)
  else
    melano <- subset(seurat_obj, Pcsk2 >= threshold)
  
  outfile <- paste0("rds_outs/", seurat_obj$orig.ident[1], "_melanotrophs.rds")
  saveRDS(melano, outfile)
  
  return(melano)
  })

corticotrophs <- pbsapply(seurat_objects_pomc, function(seurat_obj)
  {
  # Pax7 is not always detected
  if ("Pax7" %in% rownames(seurat_obj))
    cortico <- subset(seurat_obj, Pcsk2 < threshold & Pax7 < threshold)
  else
    cortico <- subset(seurat_obj, Pcsk2 < threshold)
  
  outfile <- paste0("rds_outs/", seurat_obj$orig.ident[1], "_corticotrophs.rds")
  saveRDS(cortico, outfile)

  return(cortico)
})

p <- lapply(melanotrophs, function(s){
  p <- FeaturePlot(s, "Pcsk2") +
    scale_color_viridis_c() +
    ggtitle(s$orig.ident[1])
  
  return(p)
})

do.call("grid.arrange", c(p, nrow = 3, ncol=4))

p <- lapply(corticotrophs, function(s){
  p <- FeaturePlot(s, "Pomc") +
    scale_color_viridis_c() +
    ggtitle(s$orig.ident[1])
  
  return(p)
})

do.call("grid.arrange", c(p, nrow = 3, ncol=4))

n_cells <- data.frame(dataset = sapply(seurat_objects_pomc, 
                                       function(s){s$orig.ident[1]}),
                      num_melano = sapply(melanotrophs, ncol),
                      num_cort = sapply(corticotrophs, ncol),
                      num_cells = sapply(seurat_objects, ncol))

n_cells %>% 
  mutate(perc_cort = round(num_cort/num_cells * 100, 2),
         perc_melano = round(num_melano/num_cells * 100, 2)) -> n_cells
n_cells
