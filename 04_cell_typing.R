# 04_cell_typing.R
# Clusters data, then filters them to get corticotrophs (and melanotrophs)

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(pbapply)

set.seed(12345)

plot_out_type <- "pdf" # or "pdf" or "none"
POMC_plot_folder <- "plots/Pomc_plots/"

if (!dir.exists(POMC_plot_folder)) {
  dir.create(POMC_plot_folder)
}

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

pomc_violins <- function(seurat_obj) {
  #' Plot violin plots of Pomc, Pax7, Pcsk2 and Oacyl expression
  #' @param seurat_obj: the Seurat object
  #' @return: a plot of the violin plots

  features <- c(
#    "Pomc", # Common
    "Avpr1b", "Crhr1", # Corticotrophs
    "Pax7", "Pcsk2", "Oacyl" # Melanotrophs
  ) 
  expr <- GetAssayData(seurat_obj)

  for (f in features) {
    if (!f %in% rownames(expr)) {
      expr <- rbind(expr, rep(0, ncol(expr)))
      rownames(expr)[nrow(expr)] <- f
    }
  }

  expr <- data.frame(t(as.matrix(expr[features, ]))) %>%
    mutate(cluster = Idents(seurat_obj))

  vln_plot <- expr %>%
    pivot_longer(
      cols = features,
      names_to = "feature", values_to = "expr"
    ) %>%
    ggplot(aes(x = cluster, y = expr, fill = factor(feature, levels = features))) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_fill_manual(values = c(
      # "Pomc" = "#c3e600",
      "Avpr1b" = "#c904cf",
      "Crhr1" = "#c904cf",
      "Pax7" = "#00bfc4",
      "Pcsk2" = "#00bfc4",
      "Oacyl" = "#00bfc4"
    )) +
    xlab("Cluster") +
    ylab("Expression") +
    ggtitle(paste(
      seurat_obj$author[1], seurat_obj$year[1], "-",
      seurat_obj$sex[1]
    )) +
    guides(fill = guide_legend(title = NULL)) +
    theme_bw()

  return(vln_plot)
}

get_pomc_cells <- function(seurat_obj, resolution = 0.6,
                           pomc_clusters = NULL, reduce_plot = 0.5,
                           save_rds = TRUE) {
  #' Get the cells that are later classified as corticotrophs/melanotrophs
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
  cairo_pdf("plots/PCA_elbow_plots.pdf", width = 15, height = 8)
}

do.call("grid.arrange", c(plots, nrow = 4))

if (plot_out_type != "none") {
  dev.off()
}

violins <- lapply(seurat_objects_pomc, pomc_violins)
do.call("grid.arrange", c(violins, ncol = 4))

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
  0.5 # Vennekens 2021
)

pomc_clusters <- list(
  c(4, 7, 12), # Allensworth 2021
  c(6, 11), # Cheung 2018
  3, # Fletcher 2019 F
  c(3, 10), # Fletcher 2019 M
  c(5, 6, 10), c(4, 6, 10), # Ho 2020 F/M - very noisy
  c(6, 20, 22), # Kucka 2021
  c(7, 10), # Lopez 2021
  c(1, 5), # Mayran 2019
  c(3, 9, 12, 18), c(3, 6, 15, 18), # Ruf Zamojski 2019 F/M
  4 # Vennekens 2021
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

# Violin plots of corticotrophs and melanotrophs markers

if (plot_out_type == "png") {
  png(paste(POMC_plot_folder, "Pomc_violins.png", sep = "/"), width = 12, height = 8, units = "in", res = 300)
} else if (plot_out_type == "pdf") {
  cairo_pdf(paste(POMC_plot_folder, "Pomc_violins.pdf", sep = "/"), width = 12, height = 8)
}

pomc_violins <- lapply(seurat_objects_pomc, pomc_violins)
do.call("grid.arrange", c(pomc_violins, ncol = 4))

if (plot_out_type != "none") {
  dev.off()
}

# Manually picked clusters for corticotrophs and melanotrophs, based on the violin plots
cortico_clusters <- list(
  c(7), # Allensworth 2021
  c(11), # Cheung 2018
  3, # Fletcher 2019 F
  c(3), # Fletcher 2019 M
  c(5), c(6, 10), # Ho 2020 F/M - very noisy
  c(6, 20, 22), # Kucka 2021
  c(10), # Lopez 2021
  c(1), # Mayran 2019
  c(3), c(3, 6, 15), # Ruf Zamojski 2019 F/M
  4 # Vennekens 2021
)

corticotrophs <- sapply(seq_along(seurat_objects_pomc), function(i) {
  seurat_obj <- seurat_objects_pomc[[i]]
  print(paste("Processing", seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1]))
  subset(seurat_obj, idents = cortico_clusters[[i]])
})

for (i in seq_along(corticotrophs)) {
  outfile <- paste0(
    "rds_outs/",
    seurat_objects_pomc[[i]]$orig.ident[1], "_corticotrophs.rds"
  )
  saveRDS(corticotrophs[[i]], outfile)
}

melanotrophs <- sapply(seq_along(seurat_objects_pomc), function(i) {
  seurat_obj <- seurat_objects_pomc[[i]]
  print(paste("Processing", seurat_obj$author[1], seurat_obj$year[1], "-", seurat_obj$sex[1]))

  # Some datasets have no melanotrophs, so we need to check if the setdiff is empty
  if (length(setdiff(unique(Idents(seurat_obj)), cortico_clusters[[i]])) == 0) {
    return(NULL)
  }

  return (subset(seurat_obj, idents = setdiff(unique(Idents(seurat_obj)), cortico_clusters[[i]])))
})

for (i in seq_along(melanotrophs)) {
  outfile <- paste0(
    "rds_outs/",
    seurat_objects_pomc[[i]]$orig.ident[1], "_melanotrophs.rds"
  )
  saveRDS(melanotrophs[[i]], outfile)
}


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

n_cells <- data.frame(
  dataset = sapply(
    seurat_objects_pomc,
    function(s) {
      s$orig.ident[1]
    }
  ),
  num_melano = unlist(lapply(sapply(melanotrophs, ncol), function(x) ifelse(is.null(x), 0, x))),
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
  cairo_pdf("plots/num_cells.pdf", width = 8, height = 8)
}

n_cells <- n_cells %>%
  pivot_longer(
    cols = c("perc_cort", "perc_melano"),
    names_to = "cell_type", values_to = "perc"
  ) %>%
  mutate(cell_type = factor(cell_type)) %>%
  mutate(cell_type = ifelse(cell_type == "perc_cort",
    "Corticotrophs", "Melanotrophs"
  )) %>%
  # Remove lines with 0 melanotrophs, as those correspond to datasets that
  # didn't have an intermediate zone.
  filter(!(cell_type == "Melanotrophs" & perc == 0))

n_cells %>%
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

model <- lm(perc ~ sex * cell_type, n_cells)
summary(model)
