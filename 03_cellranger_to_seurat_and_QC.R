# Imports data into Seurat and performs quality controls

library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra) # To arrange plots
library(here)

# Read datasets paths and metadata
datasets <- read.csv("datasets.csv")

load_data <- function(metadata, study) {
  matrix10x_list <- apply(metadata, 1, function(row){
      input_cellranger_dir <- row["cr_output"]
      print(paste0("Loading ", row["source"], " (", study$study_id, ") from ", 
                   input_cellranger_dir))

      res <- Read10X(paste0(getwd(), "/", input_cellranger_dir))
      })

  print("Merging matrices...")
  matrix10x <- do.call("cbind", matrix10x_list)

  # Load data into Seurat
  seurat_obj <- CreateSeuratObject(matrix10x, project = study$study_id)

  # Add metadata
  seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj,
    pattern = "^mt-"
  )

  seurat_obj[["species"]] <- metadata$species[1]
  seurat_obj[["strain"]] <- metadata$strain[1]
  seurat_obj[["sex"]] <- metadata$sex[1]
  seurat_obj[["stage"]] <- metadata$stage[1]
  seurat_obj[["age_wk"]] <- metadata$age_wk[1]
  seurat_obj[["data_source"]] <- metadata$source[1]
  
  # showWarnings = FALSE avoids warnings 
  # if the directory has already been created
  dir.create("rds_outs", showWarnings = FALSE)
  filename <- paste0("rds_outs/", study$study_id, "_raw_counts.rds")
  print(paste("Saving to", filename))
  saveRDS(object = seurat_obj, file = filename)
  
  return(seurat_obj)
}

plot_qc <- function(seurat_object,
                    min_pct_mito = NA, max_pct_mito = NA,
                    min_features = NA, max_features = NA,
                    min_counts = NA, max_counts = NA,
                    main_title = "QC plots",
                    save_to_file = FALSE) {
  #' Plots QC plots for the data
  #' This creates a 3x3 matrix of plots including
  #' Histograms of counts/cell, genes/cell, % mitochondrial genes
  #' Rank plots of the above
  #' Scatter plots of the combinations of these features
  #' @param seurat_object: The Seurat object
  #' @param min_pct_mito: Minimum percentage of mitochondrial genes
  #' @param max_pct_mito: Maximum percentage of mitochondrial genes
  #' @param min_features: Minimum number of features
  #' @param max_features: Maximum number of features
  #' @param min_counts: Minimum number of counts
  #' @param max_counts: Maximum number of counts
  #' @param main_title: Title of the plot
  #' @param save_to_file: Whether to save the plot to file (default: FALSE)

  qc_data <- data.frame(
    percent_mt = seurat_object[["percent_mt"]],
    nCount_RNA = seurat_object[["nCount_RNA"]],
    nFeature_RNA = seurat_object[["nFeature_RNA"]]
  )

  # Set limits to data limits, if not specified
  if (is.na(min_pct_mito)) {
    min_pct_mito <- min(qc_data$percent_mt)
  }
  if (is.na(max_pct_mito)) {
    max_pct_mito <- max(qc_data$percent_mt)
  }
  if (is.na(min_features)) {
    min_features <- min(qc_data$nFeature_RNA)
  }
  if (is.na(max_features)) {
    max_features <- max(qc_data$nFeature_RNA)
  }
  if (is.na(min_counts)) {
    min_counts <- min(qc_data$nCount_RNA)
  }
  if (is.na(max_counts)) {
    max_counts <- max(qc_data$nCount_RNA)
  }

  # Common theme for our plots
  if (save_to_file)
    {
    plot_theme <- theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 21))
    }  
  else
    {
    plot_theme <- theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11))
    }

  # Calculate how many points we'll discard
  qc_data %>%
    mutate(rejected = nCount_RNA >= min_counts &
      nCount_RNA <= max_counts) %>%
    select(rejected) %>%
    sum(na.rm = TRUE) -> n_reject

  perc_reject <- format(n_reject / nrow(qc_data) * 100, digits = 2, nsmall = 2)
  
  if (save_to_file)
    {
    # showWarnings = FALSE avoids warnings 
    # if the directory has already been created
    dir.create("QC_plots", showWarnings = FALSE)
    png(paste0("QC_plots/", main_title, ".png"),
        width = 1500, height = 1500, 
        pointsize = 25,
        title = main_title)
    }

  #### Histograms ####

  p1 <- ggplot(qc_data, aes(x = nCount_RNA)) +
    geom_histogram(
      binwidth = 500) +
    scale_y_log10() +
    xlab("Reads/cell") +
    ylab("Count") +
    annotate(
      geom = "text",
      x = Inf, y = Inf, # This gets mapped to the top-right corner
      label = paste(perc_reject, "% cells excluded"),
      hjust = 1.05, vjust = 1.5
    ) +
    plot_theme

  # Calculate how many points we'll discard
  qc_data %>%
    mutate(rejected = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features) %>%
    select(rejected) %>%
    sum(na.rm = TRUE) -> n_reject

  perc_reject <- format(n_reject / nrow(qc_data) * 100, digits = 2, nsmall = 2)

  p2 <- ggplot(qc_data, aes(x = nFeature_RNA)) +
    geom_histogram(
      binwidth = 100) +
    scale_y_log10() +
    xlab("Genes/cell") +
    ylab("Count") +
    annotate(
      geom = "text",
      x = Inf, y = Inf, # This gets mapped to the top-right corner
      label = paste(perc_reject, "% cells excluded"),
      hjust = 1.05, vjust = 1.5
    ) +
    plot_theme

  qc_data %>%
    mutate(rejected = percent_mt >= min_pct_mito &
      percent_mt <= max_pct_mito) %>%
    select(rejected) %>%
    sum(na.rm = TRUE) -> n_reject

  perc_reject <- format(n_reject / nrow(qc_data) * 100, digits = 2, nsmall = 2)

  p3 <- ggplot(qc_data, aes(x = percent_mt)) +
    geom_histogram(
      binwidth = 1) +
    scale_y_log10() +
    xlab("% mitochondrial genes") +
    ylab("Count") +
    annotate(
      geom = "text",
      x = Inf, y = Inf, # This gets mapped to the top-right corner
      label = paste(perc_reject, "% cells excluded"),
      hjust = 1.05, vjust = 1.5
    ) +
    plot_theme

  #### Rank plots ####
  
  # Get ranks for each cell with respect to the QC metrics
  # Note the minus because we want 1=highest and not lowest
  qc_data %>%
    mutate(CB_rank_features = rank(-nFeature_RNA)) %>%
    mutate(CB_rank_counts = rank(-nCount_RNA)) %>%
    mutate(CB_rank_pct_mt = rank(-percent_mt)) -> qc_data

  p4 <- ggplot(qc_data, aes(x = CB_rank_counts, y = nCount_RNA)) +
    geom_line() +
    geom_hline(
      yintercept = c(min_counts, max_counts),
      lty = "dotted", color = "red"
    ) +
    scale_y_log10() +
    ylab("Reads/cell") +
    xlab("Cell barcode rank") +
    plot_theme

  p5 <- ggplot(qc_data, aes(x = CB_rank_features, y = nFeature_RNA)) +
    geom_line() +
    geom_hline(
      yintercept = c(min_features, max_features),
      lty = "dotted", color = "red"
    ) +
    scale_y_log10() +
    ylab("Genes/cell") +
    xlab("Cell barcode rank") +
    plot_theme

  p6 <- ggplot(qc_data, aes(x = CB_rank_pct_mt, y = percent_mt)) +
    geom_line() +
    geom_hline(
      yintercept = c(min_pct_mito, max_pct_mito),
      lty = "dotted", color = "red"
    ) +
    scale_y_log10() +
    ylab("% mitochondrial genes") +
    xlab("Cell barcode rank") +
    plot_theme
  
  #### Scatter plots ####
  
  pt_size <- 0.1
  
  p7 <- ggplot(qc_data, aes(x = percent_mt, y = nCount_RNA)) +
    geom_point(size = pt_size) +
    geom_vline(xintercept = c(min_pct_mito, max_pct_mito), 
               col = "red",
               lty = "dotted") +
    geom_hline(yintercept = c(min_counts, max_counts), 
               col = "red",
               lty = "dotted") +
    xlab("% mitochondrial genes") +
    ylab("Reads/cell") +
    plot_theme
  
  p8 <- ggplot(qc_data, aes(x = percent_mt, y = nFeature_RNA)) +
    geom_point(size = pt_size) +
    geom_vline(xintercept = c(min_pct_mito, max_pct_mito), 
               col = "red",
               lty = "dotted") +
    geom_hline(yintercept = c(min_features, max_features), 
               col = "red",
               lty = "dotted") +
      xlab("% mitochondrial genes") +
    ylab("Genes/cell") +
    plot_theme
  
  p9 <- ggplot(qc_data, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(size = pt_size, col = rgb(0, 0, 0, 0.3)) +
    geom_vline(xintercept = c(min_counts, max_counts), 
               col = "red",
               lty = "dotted") +
    geom_hline(yintercept = c(min_features, max_features), 
               col = "red",
               lty = "dotted") +
    xlab("Reads/cell") +
    ylab("Genes/cell") +
    plot_theme

  grid.arrange(p1, p2, p3, 
               p4, p5, p6,
               p7, p8, p9,
    ncol = 3, nrow = 3,
    top = main_title
  )
  
  if (save_to_file)
    dev.off()
}

seurat_objects <- datasets %>% 
  filter(study_id != "Ho2020M") %>% 
  group_by(study_id) %>% # Group by study (M/F have different IDs, so are separated)
  group_map(~ load_data(.x, .y)) # Apply the load_data function to each group

for (i in 1:length(seurat_objects))
   {
   plot_qc(seurat_objects[[i]],
           main_title = paste(seurat_objects[[i]]$orig.ident[1], "QC"),
           save_to_file = TRUE)
    }


# Optional filtering on genes/cell

# What % of cells would we discard by further filtering on genes/cell?
# 
# nfeat <- lapply(seurat_objects, function(x){
#               data.frame(nFeatures = x$nFeature_RNA, 
#                          dataset = x$orig.ident[1])
#         })
# nfeat <- do.call("rbind", nfeat)
# 
# nfeat %>% 
#   group_by(dataset) %>% 
#   summarise(num_cells = length(nFeatures),
#             perc_less_50 = sum(nFeatures < 50)/sum(nFeatures>0)*100,
#             perc_less_100 = sum(nFeatures < 100)/sum(nFeatures>0)*100,
#             perc_less_200 = sum(nFeatures < 200)/sum(nFeatures>0)*100,
#             perc_less_500 = sum(nFeatures < 500)/sum(nFeatures>0)*100)
# 
# seurat_objects <- lapply(seurat_objects, function(x)
#   {
#   x %>% 
#     subset(nFeature_RNA > 100)
#   })
# 
# # Re-plot QC
# for (i in 1:length(seurat_objects))
#   {
#   plot_qc(seurat_objects[[i]],
#           main_title = paste(seurat_objects[[i]]$orig.ident[1], "QC"),
#           save_to_file = TRUE)
#   }