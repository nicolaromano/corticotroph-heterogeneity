# Imports data into Seurat and performs quality controls

library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pbapply)

# Read datasets paths and metadata
datasets <- read.csv("datasets.csv")

load_data <- function(metadata, study) {
  #' Loads data from CellRanger output, adds metadata, filters and saves
  #' the resulting Seurat object in an RDS file in the rds_outs folder
  #' @param metadata: the metadata for the datasets
  #' @param study: the ID of the study we are loading
  #' Parameters for filtering. Cells are filtered if they have less than min_quantile_xxx or min_xxx, whichever is higher.
  #' Similarly for max_quantile_xxx and max_xxx.
  #' @return: a Seurat object with the data and metadata

  matrix10x_list <- apply(metadata, 1, function(row) {
    input_cellranger_dir <- row["cr_output"]
    print(paste0(
      "Loading ", row["source"], " (", study$study_id, ") from ",
      input_cellranger_dir
    ))

    res <- Read10X(paste0(getwd(), "/", input_cellranger_dir))

    return(res)
  })

  print("Merging matrices...")
  matrix10x <- do.call("cbind", matrix10x_list)

  # Load data into Seurat
  seurat_obj <- CreateSeuratObject(matrix10x, project = study$study_id)

  # Add metadata
  if (metadata$species[1] == "Mouse") {
    seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj,
      pattern = "^mt-"
    )
  } else {
    seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj,
      pattern = "^Mt-"
    )
  }
  seurat_obj[["species"]] <- metadata$species[1]
  seurat_obj[["strain"]] <- metadata$strain[1]
  seurat_obj[["sex"]] <- metadata$sex[1]
  seurat_obj[["stage"]] <- metadata$stage[1]
  seurat_obj[["age_wk"]] <- metadata$age_wk[1]
  seurat_obj[["data_source"]] <- metadata$source[1]

  # Filtering cells

  # We use a relatively loose filtering here.
  # Filter anything greater than 75th percentile + 3 * SD for counts and features
  # Filter anything with <100 genes or <300 counts
  # and anything with >80% mitochondrial genes
  low_counts <- 300
  high_counts <- quantile(seurat_obj$nCount_RNA, 0.75) + 3 * sd(seurat_obj$nCount_RNA)
  low_features <- 100
  high_features <- quantile(seurat_obj$nFeature_RNA, 0.75) + 3 * sd(seurat_obj$nFeature_RNA)
  low_mt <- 0
  high_mt <- 80

  print(paste("Loaded", ncol(seurat_obj), "cells"))

  # Set the limits for the plots so that the post-filtering plots are comparable to the pre-filtering ones
  lim_counts <- c(0, max(seurat_obj$nCount_RNA) * 1.1)
  lim_features <- c(0, max(seurat_obj$nFeature_RNA) * 1.1)
  lim_mt <- c(0, max(seurat_obj$percent_mt) * 1.1)

  plot_qc(seurat_obj,
    main_title = paste(seurat_obj$orig.ident[1], "QC_pre_filtering"),
    save_to_file = TRUE,
    lim_counts = lim_counts,
    lim_features = lim_features,
    lim_mt = lim_mt
  ) 

  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= low_features & nFeature_RNA <= high_features &
      nCount_RNA >= low_counts & nCount_RNA <= high_counts &
      percent_mt >= low_mt & percent_mt <= high_mt
  )

  print(paste("Filtered to", ncol(seurat_obj), "cells"))

  plot_qc(seurat_obj,
    main_title = paste(seurat_obj$orig.ident[1], "QC_post_filtering"),
    save_to_file = TRUE,
    lim_counts = lim_counts,
    lim_features = lim_features,
    lim_mt = lim_mt
  )

  # showWarnings = FALSE avoids warnings
  # if the directory has already been created
  dir.create("rds_outs", showWarnings = FALSE)
  filename <- paste0("rds_outs/", study$study_id, "_raw_counts.rds")
  print(paste("Saving to", filename))
  saveRDS(object = seurat_obj, file = filename)

  return(seurat_obj)
}

plot_qc <- function(seurat_object,
                    main_title = "QC plots",
                    save_to_file = FALSE,
                    lim_features = NULL,
                    lim_counts = NULL,
                    lim_mt = NULL) {                    
  #' Plots QC plots for the data
  #' This creates a 3x3 matrix of plots including
  #' Histograms of counts/cell, genes/cell, % mitochondrial genes
  #' Rank plots of the above
  #' Scatter plots of the combinations of these features
  #' @param seurat_object: The Seurat object
  #' @param main_title: Title of the plot
  #' @param save_to_file: Whether to save the plot to file (default: FALSE)
  #' @param lim_features: axis limits for the genes/cell plot. If NULL, the limits are automatically calculated
  #' @param lim_counts: axis limits for the counts/cell plot. If NULL, the limits are automatically calculated  
  #' @param lim_mt: axis limits for the % mitochondrial genes plot. If NULL, the limits are automatically calculated

  qc_data <- data.frame(
    percent_mt = seurat_object[["percent_mt"]],
    nCount_RNA = seurat_object[["nCount_RNA"]],
    nFeature_RNA = seurat_object[["nFeature_RNA"]]
  )

  # Common theme for our plots
  if (save_to_file) {
    plot_theme <- theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 21)
    )
  } else {
    plot_theme <- theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
  }

  if (save_to_file) {
    # showWarnings = FALSE avoids warnings
    # if the directory has already been created
    dir.create("plots/QC_plots", showWarnings = FALSE)
    png(paste0("plots/QC_plots/", main_title, ".png"),
      width = 1500, height = 1500,
      pointsize = 25,
      title = main_title
    )
  }

  #### Histograms ####

  p1 <- ggplot(qc_data, aes(x = nCount_RNA)) +
    geom_histogram(
      binwidth = 500
    ) +
    scale_y_log10() +
    xlab("Reads/cell") +
    ylab("Count") +
    plot_theme

  if (!is.null(lim_counts)) {
    p1 <- p1 + xlim(lim_counts)
  }

  p2 <- ggplot(qc_data, aes(x = nFeature_RNA)) +
    geom_histogram(
      binwidth = 100
    ) +
    scale_y_log10() +
    xlab("Genes/cell") +
    ylab("Count") +
    plot_theme

  if (!is.null(lim_features)) {
    p2 <- p2 + xlim(lim_features)
  }

  p3 <- ggplot(qc_data, aes(x = percent_mt)) +
    geom_histogram(
      binwidth = 1
    ) +
    scale_y_log10() +
    xlab("% mitochondrial genes") +
    ylab("Count") +
    plot_theme

  if (!is.null(lim_mt)) {
    p3 <- p3 + xlim(lim_mt)
  }

  #### Rank plots ####

  # Get ranks for each cell with respect to the QC metrics
  # Note the minus because we want 1=highest and not lowest
  qc_data %>%
    mutate(CB_rank_features = rank(-nFeature_RNA)) %>%
    mutate(CB_rank_counts = rank(-nCount_RNA)) %>%
    mutate(CB_rank_pct_mt = rank(-percent_mt)) -> qc_data

  p4 <- ggplot(qc_data, aes(x = CB_rank_counts, y = nCount_RNA)) +
    geom_line() +
    scale_y_log10() +
    ylab("Reads/cell") +
    xlab("Cell barcode rank") +
    plot_theme

  p5 <- ggplot(qc_data, aes(x = CB_rank_features, y = nFeature_RNA)) +
    geom_line() +
    scale_y_log10() +
    ylab("Genes/cell") +
    xlab("Cell barcode rank") +
    plot_theme

  p6 <- ggplot(qc_data, aes(x = CB_rank_pct_mt, y = percent_mt)) +
    geom_line() +
    scale_y_log10() +
    ylab("% mitochondrial genes") +
    xlab("Cell barcode rank") +
    plot_theme

  #### Scatter plots ####

  pt_size <- 0.1

  p7 <- ggplot(qc_data, aes(x = percent_mt, y = nCount_RNA)) +
    geom_point(size = pt_size) +
    xlab("% mitochondrial genes") +
    ylab("Reads/cell") +
    plot_theme
  
  if (!is.null(lim_mt)) {
    p7 <- p7 + xlim(lim_mt)
  }
  if (!is.null(lim_counts)) {
    p7 <- p7 + ylim(lim_counts)
  }

  p8 <- ggplot(qc_data, aes(x = percent_mt, y = nFeature_RNA)) +
    geom_point(size = pt_size) +
    xlab("% mitochondrial genes") +
    ylab("Genes/cell") +
    plot_theme

  if (!is.null(lim_mt)) {
    p8 <- p8 + xlim(lim_mt)
  }
  if (!is.null(lim_features)) {
    p8 <- p8 + ylim(lim_features)
  }
  
  p9 <- ggplot(qc_data, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(size = pt_size, col = rgb(0, 0, 0, 0.3)) +
    xlab("Reads/cell") +
    ylab("Genes/cell") +
    plot_theme

  if (!is.null(lim_counts)) {
    p9 <- p9 + xlim(lim_counts)
  }
  if (!is.null(lim_features)) {
    p9 <- p9 + ylim(lim_features)
  }

  grid.arrange(p1, p2, p3,
    p4, p5, p6,
    p7, p8, p9,
    ncol = 3, nrow = 3,
    top = main_title
  )

  if (save_to_file) {
    dev.off()
  }
}

seurat_objects <- datasets %>%
  # Group by study (M/F have different IDs, so are separated)
  group_by(study_id) %>%
  # Apply the load_data function to each group
  # Filter with the default settings (see load_data function)
  group_map(~ load_data(.x, .y))

# Optionally, load the raw counts directly from the RDS files
# filenames <- list.files("rds_outs", pattern = "raw_counts.rds", full.names = TRUE)
# seurat_objects <- filenames %>%
#   pblapply(readRDS)

# Print number of cells
num_cells <- data.frame(
  dataset = sapply(seurat_objects, function(s) {
    as.character(s$orig.ident[1])
  }),
  cells = sapply(seurat_objects, ncol)
)

num_cells

# Scale with SCT -> PCA -> UMAP
seurat_objects_SCT <- sapply(seurat_objects, function(s) {
  s %>%
    SCTransform() %>%
    FindVariableFeatures() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:20) -> seuratobj

  outfile <- paste0("rds_outs/", s$orig.ident[1], "_SCT.rds")
  saveRDS(seuratobj, file = outfile)

  return(seuratobj)
})
