# Imports data into Seurat and performs quality controls

library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra) # To arrange plots
library(here)

# Read datasets paths and metadata
datasets <- read.csv("datasets.csv")

load_data <- function(metadata) {
  input_cellranger_dir <- metadata$cr_output
  project_name <- metadata$id

  print(paste("Loading", project_name, "from", input_cellranger_dir))
  # Load data into Seurat
  seurat_obj <- CreateSeuratObject(Read10X(paste0(
    getwd(), "/",
    input_cellranger_dir
  )),
  project = project_name
  )

  # Add metadata
  seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj,
    pattern = "^mt-"
  )

  seurat_obj[["species"]] <- metadata$species
  seurat_obj[["strain"]] <- metadata$strain
  seurat_obj[["sex"]] <- metadata$sex
  seurat_obj[["stage"]] <- metadata$stage
  seurat_obj[["age_wk"]] <- metadata$age_wk

  return(seurat_obj)
}

plot_qc <- function(seurat_object, save_to_file = FALSE,
                    min_pct_mito = NA, max_pct_mito = NA,
                    min_features = NA, max_features = NA,
                    min_counts = NA, max_counts = NA) {
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
  plot_theme <- theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

  # Calculate how many points we'll discard
  qc_data %>%
    mutate(rejected = nCount_RNA >= min_counts &
      nCount_RNA <= max_counts) %>%
    select(rejected) %>%
    sum(na.rm = TRUE) -> n_reject

  perc_reject <- format(n_reject / nrow(qc_data) * 100, digits = 2, nsmall = 2)
  
  #### Histograms ####

  p1 <- ggplot(qc_data, aes(x = nCount_RNA)) +
    geom_histogram(
      binwidth = 500,
      aes(fill = nCount_RNA >= min_counts &
        nCount_RNA <= max_counts)
    ) +
    scale_fill_manual(
      values = c(rgb(0.8, 0.1, 0.1), rgb(0.3, 0.3, 0.3)),
      guide = "none"
    ) +
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
      binwidth = 100,
      aes(fill = nFeature_RNA >= min_features &
        nFeature_RNA <= max_features)
    ) +
    scale_fill_manual(
      values = c(rgb(0.8, 0.1, 0.1), rgb(0.3, 0.3, 0.3)),
      guide = "none"
    ) +
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
      binwidth = 1,
      aes(fill = percent_mt >= min_pct_mito &
        percent_mt <= max_pct_mito)
    ) +
    scale_fill_manual(
      values = c(rgb(0.8, 0.1, 0.1), rgb(0.3, 0.3, 0.3)),
      guide = "none"
    ) +
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
    ncol = 3, nrow = 3
  )
}

# seurat_object <- load_data(datasets[3,])
plot_qc(seurat_object,
  min_counts = 500,
  min_features = 200
)

#
# # tiff(filename = paste0("QC plots", project_name, "pre-filtering", sep=" - "),
# # width = 800, height = 800, compression = "lzwm")
#
# ggarrange(p1, p2, p3, ncol = 2, nrow = 2) +
#   ggtitle(paste("QC plots", project_name, "pre-filtering", sep = " - "))
# # dev.off()
