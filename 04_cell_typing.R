# 04_cell_typing.R
# Clusters data, then filters them to get corticotrophs (and melanotrophs)

library(Seurat)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(pbapply)

set.seed(12345)

plot_out_type <- "none" # or "pdf" or "none"

#### Various helper functions ####
show_markers <- function(seurat_obj, hormonal = TRUE) {
  #' Plots [non]-hormonal features for the given dataset
  #' @param seurat_obj Seurat object
  #' @param hormonal Whether to plot features of hormonal cell
  #' types or those of other cell types
  #' @return A plot of the features

  if (hormonal) {
    features <- c(
      "Gh", # Somatotrophs
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

  # Note combine = FALSE to return a list of plots
  p <- FeaturePlot(seurat_obj, features, pt.size = 1, combine = FALSE)

  p <- lapply(p, function(x) {
    x +
      scale_color_viridis_c() +
      xlab(expression(UMAP[1])) +
      ylab(expression(UMAP[2]))
  })

  do.call(
    "grid.arrange",
    c(p,
      ncol = 3,
      top = paste(seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1])
    )
  )
}

plot_features_histo <- function(seurat_obj, feature, find_threshold = TRUE,
                                log_y_axis = TRUE) {
  if (!feature %in% rownames(seurat_obj)) # Feature does not exist
    {
      p <- ggplot() +
        ggtitle(paste(
          seurat_obj$author[1], seurat_obj$year[1], "-",
          seurat_obj$sex[1]
        ), subtitle = "Not found")

      return(p)
    }

  expr <- data.frame(expr = GetAssayData(seurat_obj)[feature, ])

  p <- ggplot(expr, aes(x = expr)) +
    geom_histogram(binwidth = 0.1) +
    xlim(0, 10) +
    xlab(paste(feature, "expression")) +
    ggtitle(paste(
      seurat_obj$author[1],
      seurat_obj$year[1], "-", seurat_obj$sex[1]
    ))

  if (log_y_axis) {
    p <- p + scale_y_log10()
  }

  if (find_threshold) {
    thr <- otsu_thresh(expr$expr, levels = 500)
    p <- p + geom_vline(xintercept = thr, col = "red", lty = "dashed")
  }

  p
}

otsu_thresh <- function(values, levels = 100) {
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
  var <- w1 * w2 * (m2 / w2 - m1 / w1)^2
  maxi <- which(var == max(var, na.rm = TRUE))

  return((mids[maxi[1]] + mids[maxi[length(maxi)]]) / 2)
}

type_pomc_cells <- function(seurat_obj) {
  #' Type cells based on a Otsu's threshold on the expression of a feature
  #' @param seurat_obj: the Seurat object
  #' @return: a plot of the cell types

  expr <- GetAssayData(seurat_obj)
  thr_Pomc <- otsu_thresh(expr["Pomc", ])
  if ("Pax7" %in% rownames(expr))
    thr_Pax7 <- otsu_thresh(expr["Pax7", ])
  if ("Pcsk2" %in% rownames(expr))
    thr_Pcsk2 <- otsu_thresh(expr["Pcsk2", ])

  cells_bc <- Cells(subset(seurat_obj, Pomc > thr_Pomc))
  print(paste0(
    round(length(cells_bc) / length(Cells(seurat_obj)) * 100, 2),
    "% cells typed as corticotrophs/melanotrophs"
  ))

  Pomc_pos <- expr["Pomc", ] > thr_Pomc
  if ("Pax7" %in% rownames(expr) & "Pcsk2" %in% rownames(expr)) {
    Pax7_Pcsk2_pos <- expr["Pax7", ] > thr_Pax7 | expr["Pcsk2", ] > thr_Pcsk2
  } else if ("Pax7" %in% rownames(expr)) {
    Pax7_Pcsk2_pos <- expr["Pax7", ] > thr_Pax7
  } else if ("Pcsk2" %in% rownames(expr)) {
    Pax7_Pcsk2_pos <- expr["Pcsk2", ] > thr_Pcsk2
  } else {
    Pax7_Pcsk2_pos <- rep(0, length(cells_bc))
  }

  seurat_obj$finalIdent <- ifelse(Pomc_pos + Pax7_Pcsk2_pos == 0, "Other",
    ifelse(Pomc_pos + Pax7_Pcsk2_pos == 1, "Corticotroph",
      "Melanotroph"
    )
  )
  
  p <- FeaturePlot(seurat_obj, "finalIdent", pt.size = 0.2) +
    scale_color_manual(name = "Cell type", values=c("#c3e600", "#c904cf", "lightgray")) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    theme(legend.position = "none") +
    ggtitle(paste(seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1]))
    
  return(p)
}

get_pomc_cells <- function(seurat_obj, resolution = 0.6,
                           pomc_clusters = NULL, reduce_plot = 0.5,
                           save_rds = TRUE) {
  #' Get the cells that are classified as corticotrophs/melanotrophs
  #' @param seurat_obj: the Seurat object
  #' @param resolution: the resolution to use for clustering
  #' @param pomc_clusters: the clusters to consider as corticotrophs/melanotrophs
  #' @param reduce_plot: the fraction of cells to plot (to speed up plotting)
  #' @param save_rds: whether to save the Seurat object as an RDS file
  #' @return: a Seurat object with only the cells classified as corticotrophs/melanotrophs

  seurat_obj %>%
    FindNeighbors(k.param = 20) %>%
    FindClusters(resolution = resolution) -> seurat_obj

  to_plot <- seurat_obj[, sample(colnames(seurat_obj),
    size = round(reduce_plot * ncol(seurat_obj))
  )]

  p1 <- DimPlot(to_plot, label = TRUE, pt.size = .2) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(seurat_obj$orig.ident[1]) +
    theme(legend.position = "none")

  p2 <- FeaturePlot(to_plot, feature = "Pomc", pt.size = .2) +
    scale_color_viridis_c() +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2]))

  p3 <- VlnPlot(seurat_obj, features = "Pomc") + NoLegend()

  if (plot_out_type == "png") {
    png(paste0("plots/Pomc_plots/", seurat_obj$orig.ident[1], "_all_pit_Pomc.png"),
      width = 10, height = 10, units = "in", res = 300
    )
  } else if (plot_out_type == "pdf") {
    pdf(paste0("plots/Pomc_plots/", seurat_obj$orig.ident[1], "_all_pit_Pomc.pdf"),
      width = 10, height = 10
    )
  }
  # Arrange with p1 and p2 on top, side by side, and p3 on the bottom
  grid.arrange(p1, p2, p3, ncol = 2, layout_matrix = rbind(c(1, 2), c(3, 3)))

  if (plot_out_type != "none") {
    dev.off()
  }

  if (save_rds) {
    seurat_obj %>%
      subset(idents = pomc_clusters) %>%
      RunUMAP(dims = 1:20) -> seurat_obj_pomc

    saveRDS(seurat_obj_pomc, file = paste0(
      "rds_outs/",
      seurat_obj_pomc$orig.ident[1],
      "_pomc.rds"
    ))
    return(seurat_obj_pomc)
  }

  return(seurat_obj)
}

datasets <- read.csv("datasets.csv")

filenames <- dir("rds_outs", pattern = "SCT.rds", full.names = TRUE)

# Read all datasets
seurat_objects <- pblapply(filenames, readRDS)

names(seurat_objects) <- sapply(seurat_objects, function(s) {
  s$orig.ident[1]
})

# We can use this to explore the datasets
# show_markers(seurat_objects[[2]])

# PCA elbow plots
plots <- lapply(seurat_objects, function(s) {
  p <- ElbowPlot(s, 30) +
    ggtitle(paste(s$author[1], s$year[1], "-", s$sex[1]))
})

if (plot_out_type == "png") {
  png("plots/PCA_elbow_plots.png", width = 15, height = 8, units = "in", res = 300)
} else if (plot_out_type == "pdf") {
  pdf("plots/PCA_elbow_plots.pdf", width = 15, height = 8)
}

do.call("grid.arrange", c(plots, nrow = 4))

if (plot_out_type != "none") {
  dev.off()
}

#### POMC HISTOGRAM #####
pomc_hist <- lapply(seurat_objects, plot_features_histo, "Pomc",
  find_threshold = TRUE, log_y_axis = TRUE
)

if (plot_out_type == "png") {
  png("plots/pomc_hist.png", width = 15, height = 8, units = "in", res = 300)
} else if (plot_out_type == "pdf") {
  pdf("plots/pomc_hist.pdf", width = 15, height = 8)
}

do.call("grid.arrange", c(pomc_hist, ncol = 4, nrow = 3))

if (plot_out_type != "none") {
  dev.off()
}

pomc_cells <- lapply(seurat_objects, type_pomc_cells)

if (plot_out_type == "png") {
  png("plots/pomc_cells.png", width = 15, height = 10, units = "in", res = 300)
} else if (plot_out_type == "pdf") {
  pdf("plots/pomc_cells.pdf", width = 15, height = 10)
}

do.call("grid.arrange", c(pomc_cells,
  top = "POMC+ cells",
  ncol = 4
))

if (plot_out_type != "none") {
  dev.off()
}

#### CLUSTERING ####

# These are manually picked as not possible to use the same for everyone
resolutions <- c(
  0.5, # Allensworth 2021
  0.5, # Cheung 2018
  0.3, 0.3, # Fletcher 2019 F/M
  1, 1, # Ho 2020 F/M
  0.5, # Kucka 2021
  0.5, # Lopez 2021
  0.5, # Mayran 2019
  0.3, 0.3, # Ruf Zamojski 2019 F/M
  0.5
) # Vennekens 2021

pomc_clusters <- list(
  c(5, 7, 12), # Allensworth 2021
  c(6, 11), # Cheung 2018
  2, c(3, 11), # Fletcher 2019 F/M
  c(4, 5, 6), c(3, 6, 10), # Ho 2020 F/M
  c(7, 21, 22, 23), # Kucka 2021
  5, # Lopez 2021
  c(3, 6), # Mayran 2019
  c(4, 5, 18), c(1, 16, 18), # Ruf Zamojski 2019 F/M
  5 # Vennekens 2021
)

seurat_objects_pomc <- sapply(seq_along(seurat_objects), function(i) {
  seurat_obj <- seurat_objects[[i]]
  print(paste("Processing", seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1]))
  get_pomc_cells(seurat_obj,
    resolution = resolutions[i],
    pomc_clusters = pomc_clusters[[i]],
    # Set this to a small number (0.1-0.2) while determining the resolution
    # and to 1 when actually plotting the final results
    reduce_plot = 1,
    # Set this to FALSE while determining the resolution and
    # pomc_clusters, then set it to TRUE to save the RDS files
    save_rds = TRUE
  )
})

#### SPLIT POMC CELLS INTO CORTICOTROPHS AND MELANOTROPHS ####

# Load files (to start directly from here)
filenames <- dir(path = "rds_outs/", pattern = "pomc.rds")
seurat_objects_pomc <- pbsapply(paste0("rds_outs/", filenames), readRDS)

names(seurat_objects_pomc) <- sapply(seurat_objects_pomc, function(s) {
  s$orig.ident[1]
})

pomc_plots <- lapply(seurat_objects_pomc, function(seurat_obj) {
  DimPlot(seurat_obj, label = TRUE) +
    ggtitle(paste(seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1])) +
    theme(legend.position = "none")
})

do.call("grid.arrange", c(pomc_plots, ncol = 4))

melano_plots <- lapply(seurat_objects_pomc, function(seurat_obj) {
  melano_markers <- c("Pcsk2", "Pax7")
  melano_markers <- melano_markers[melano_markers %in% rownames(seurat_obj)]

  if (length(melano_markers) > 1) {
    melano_mark_expr <- data.frame(t(as.matrix(GetAssayData(seurat_obj)[melano_markers, ])))
  } else {
    melano_mark_expr <- data.frame(as.matrix(GetAssayData(seurat_obj)[melano_markers, ]))
  }
  umap_coord <- data.frame(seurat_obj@reductions$umap@cell.embeddings)
  p <- ggplot(umap_coord, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(col = rowSums(melano_mark_expr) > 0), size = 0.5) +
    scale_color_manual(values = c("lightgray", "navy"), guide = "none") +
    ggtitle(seurat_obj$orig.ident[1])
  return(p)
})

if (plot_out_type == "png") {
  png("plots/melanotrophs.png", width = 800, height = 800)
} else if (plot_out_type == "pdf") {
  pdf("plots/melanotrophs.pdf", width = 8, height = 8)
}

do.call("grid.arrange", c(melano_plots,
  nrow = 3,
  top = "Melanotrophs (Pcsk2+ or Pax7+)"
))

if (plot_out_type != "none") {
  dev.off()
}

threshold <- 0.1

melanotrophs <- pbsapply(seurat_objects_pomc, function(seurat_obj) {
  # Pax7 is not always detected
  if ("Pax7" %in% rownames(seurat_obj)) {
    melano <- subset(seurat_obj, Pcsk2 >= threshold | Pax7 >= threshold)
  } else {
    melano <- subset(seurat_obj, Pcsk2 >= threshold)
  }

  outfile <- paste0(
    "rds_outs/",
    seurat_obj$orig.ident[1], "_melanotrophs.rds"
  )
  saveRDS(melano, outfile)

  return(melano)
})

corticotrophs <- pbsapply(seurat_objects_pomc, function(seurat_obj) {
  # Pax7 is not always detected
  if ("Pax7" %in% rownames(seurat_obj)) {
    cortico <- subset(seurat_obj, Pcsk2 < threshold & Pax7 < threshold)
  } else {
    cortico <- subset(seurat_obj, Pcsk2 < threshold)
  }

  outfile <- paste0(
    "rds_outs/", seurat_obj$orig.ident[1],
    "_corticotrophs.rds"
  )
  saveRDS(cortico, outfile)

  return(cortico)
})

cort_melano_plots <- lapply(seurat_objects_pomc, function(seurat_obj) {
  markers <- c("Pomc", "Pcsk2", "Pax7")
  markers <- markers[markers %in% rownames(seurat_obj)]

  if (length(markers) > 1) {
    mark_expr <- data.frame(t(as.matrix(GetAssayData(seurat_obj)[markers, ])))
  } else {
    mark_expr <- data.frame(as.matrix(GetAssayData(seurat_obj)[markers, ]))
  }

  umap_coord <- data.frame(seurat_obj@reductions$umap@cell.embeddings)
  p <- ggplot(umap_coord, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(col = rowSums(melano_mark_expr) > 0), size = 0.5) +
    scale_color_manual(values = c("lightgray", "navy"), guide = "none") +
    ggtitle(seurat_obj$orig.ident[1])
  return(p)
})



# Can load processed files here directly, if needed
# melanotrophs <- pbsapply(dir(
#   path = "rds_outs/",
#   pattern = "_melanotrophs.rds",
#   full.names = TRUE
# ), readRDS)

# corticotrophs <- pbsapply(dir(
#   path = "rds_outs/",
#   pattern = "_corticotrophs.rds",
#   full.names = TRUE
# ), readRDS)

p <- lapply(melanotrophs, function(s) {
  p <- FeaturePlot(s, "Pcsk2") +
    scale_color_viridis_c() +
    ggtitle(s$orig.ident[1])

  return(p)
})

do.call("grid.arrange", c(p, ncol = 4))

p <- lapply(corticotrophs, function(s) {
  p <- FeaturePlot(s, "Pomc") +
    scale_color_viridis_c() +
    ggtitle(s$orig.ident[1])

  return(p)
})

do.call("grid.arrange", c(p, ncol = 4))

n_cells <- data.frame(
  dataset = sapply(
    seurat_objects_pomc,
    function(s) {
      s$orig.ident[1]
    }
  ),
  num_melano = sapply(melanotrophs, ncol),
  num_cort = sapply(corticotrophs, ncol),
  num_cells = sapply(seurat_objects, ncol)
)

n_cells %>%
  remove_rownames() %>%
  mutate(
    perc_cort = round(num_cort / num_cells * 100, 2),
    perc_melano = round(num_melano / num_cells * 100, 2),
    sex = str_sub(dataset, -1, -1)
  ) -> n_cells
n_cells

if (plot_out_type == "png") {
  png("plots/num_cells.png", width = 8, height = 8, units = "in", res = 200)
} else if (plot_out_type == "pdf") {
  pdf("plots/num_cells.pdf", width = 8, height = 8)
}

n_cells %>%
  pivot_longer(
    cols = c("perc_cort", "perc_melano"),
    names_to = "cell_type", values_to = "perc"
  ) %>%
  mutate(cell_type = factor(cell_type)) %>%
  mutate(cell_type = ifelse(cell_type == "perc_cort",
    "Corticotrophs", "Melanotrophs"
  )) %>%
  ggplot(aes(
    x = sex, y = perc, fill = cell_type
  )) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(
    "",
    values = c("white", "lightgray")
  ) +
  ylab("% cells") +
  xlab("Sex") +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 16)
  )

if (plot_out_type != "none") {
  dev.off()
}

t.test(num_cort ~ sex, n_cells)
