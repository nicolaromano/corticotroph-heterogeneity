# 03_cellranger_to_seurat_and_QC.R
# Imports data into Seurat and performs quality controls
# Inputs: CellRanger output folders
# Outputs: RDS files with Seurat objects (rds_outs/dataset_SCT.rds), QC plots (plots/QC_plots)

library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pbapply)

plot_out_type <- "pdf" # "png" or "pdf" or "none"
load_CR_matrices <- FALSE # Set to TRUE to reload data from CellRanger output

# Read datasets paths and metadata
datasets <- read.csv("datasets.csv")

load_data <- function(metadata, study) {
  #' Loads data from CellRanger output, adds metadata, filters and saves
  #' the resulting Seurat object in an RDS file in the rds_outs folder
  #' @param metadata: the metadata for the datasets
  #' @param study: the ID of the study we are loading
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
  seurat_obj[["author"]] <- metadata$author[1]
  seurat_obj[["year"]] <- metadata$year[1]

  # Filtering cells

  # We use a relatively loose filtering here.
  # Filter anything greater than 75th percentile + 3 * SD for counts and features,
  # anything below 25th percentile - 3 * SD for counts and features (or 300 and 100 respectively, whichever is higher),
  # and anything with >80% mitochondrial genes
  low_counts <- max(500, quantile(seurat_obj$nCount_RNA, 0.25) - 3 * sd(seurat_obj$nCount_RNA))
  high_counts <- quantile(seurat_obj$nCount_RNA, 0.75) + 3 * sd(seurat_obj$nCount_RNA)
  low_features <- max(300, quantile(seurat_obj$nCount_RNA, 0.25) - 3 * sd(seurat_obj$nCount_RNA))
  high_features <- quantile(seurat_obj$nFeature_RNA, 0.75) + 3 * sd(seurat_obj$nFeature_RNA)
  low_mt <- 0
  high_mt <- 70

  print("Filtering cells with:")
  print(paste("Counts:", low_counts, high_counts))
  print(paste("Features:", low_features, high_features))
  print(paste("Mitochondrial genes:", low_mt, high_mt))

  print(paste("Loaded", ncol(seurat_obj), "cells"))

  # Set the limits for the plots so that the post-filtering plots are comparable to the pre-filtering ones
  lim_counts <- c(0, max(seurat_obj$nCount_RNA) * 1.1)
  lim_features <- c(0, max(seurat_obj$nFeature_RNA) * 1.1)
  lim_mt <- c(0, max(seurat_obj$percent_mt) * 1.1)

  plot_qc(seurat_obj,
    main_title = paste(seurat_obj$author, seurat_obj$year, seurat_obj$sex, "QC pre-filtering"),
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
    main_title = paste(seurat_obj$author, seurat_obj$year, seurat_obj$sex, "QC post-filtering"),
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
    if (plot_out_type == "pdf") {
      pdf(paste0("plots/QC_plots/", main_title, ".pdf"),
        width = 15, height = 15,
        pointsize = 25,
        title = main_title
      )
    } else if (plot_out_type == "png") {
      png(paste0("plots/QC_plots/", main_title, ".png"),
        width = 1500, height = 1500,
        pointsize = 25,
        title = main_title
      )
    }
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

if (load_CR_matrices) {
  seurat_objects <- datasets %>%
    # Group by study (M/F have different IDs, so are separated)
    group_by(study_id) %>%
    # Apply the load_data function to each group
    # Filter with the default settings (see load_data function)
    group_map(~ load_data(.x, .y))
} else {
  # Load the raw counts directly from the RDS files
  filenames <- list.files("rds_outs", pattern = "raw_counts.rds", full.names = TRUE)
  seurat_objects <- filenames %>%
    pblapply(readRDS)
}

# Print number of cells
cell_num <- data.frame(
  cells = sapply(seurat_objects, function(s) {
    ncol(s)
  })
) %>%
  mutate(reported = c(18301, 13663, 3334, 3562, 1340, 1440, 15876, 9879, 9269, 21648, 10311, 7977)) %>%
  mutate(author = sapply(seurat_objects, function(s) {
    s$author[1]
  })) %>%
  mutate(sex = sapply(seurat_objects, function(s) {
    s$sex[1]
  })) %>%
  mutate(sex = factor(sex, levels = c("M", "F"))) %>%
  mutate(year = sapply(seurat_objects, function(s) {
    s$year[1]
  })) %>%
  mutate(perc_reported = round(cells / reported * 100, 2))

# Stats - is the ratio reported vs. published different from 1?
# Does sex affect the ratio?
model <- lm(I((cells / reported) - 1) ~ sex, data = cell_num)
summary(model)

if (plot_out_type == "pdf") {
  pdf("plots/cell_number_perc_reported.pdf",
    width = 10, height = 8
  )
} else if (plot_out_type == "png") {
  png("plots/cell_number_perc_reported.png",
    width = 1000, height = 500
  )
}

ggplot(
  cell_num, aes(
    x = author, y = perc_reported,
    fill = factor(sex, levels = c("F", "M"))
  )
) +
  geom_col(width = 0.8, position = position_dodge(preserve = "single")) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  ylab("Cells\n(% published num)") +
  xlab("") +
  scale_y_continuous(limits = c(0, 300), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

if (plot_out_type != "none") {
  dev.off()
}

if (plot_out_type == "pdf") {
  pdf("plots/cell_number.pdf", width = 10, height = 8)
} else if (plot_out_type == "png") {
  png("plots/cell_number.png", width = 1000, height = 500)
}

ggplot(cell_num, aes(
  x = author, y = cells / 1e3,
  fill = factor(sex, levels = c("F", "M"))
)) +
  geom_col(position = position_dodge(preserve = "single")) +
  ylim(0, max(cell_num$cells)) +
  ylab(expression("Cells (x " ~ 10^3 ~ ")")) +
  xlab("") +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

if (plot_out_type != "none") {
  dev.off()
}

# Plot QC for all studies
all_qc_data <- data.frame(
  counts = sapply(
    seurat_objects,
    function(s) {
      s$nCount_RNA
    }
  ) %>% unlist(),
  genes = sapply(seurat_objects, function(s) {
    s$nFeature_RNA
  }) %>% unlist(),
  mt_perc = sapply(seurat_objects, function(s) {
    s$percent_mt
  }) %>% unlist(),
  dataset = sapply(seurat_objects, function(s) {
    s$author
  }) %>% unlist(),
  sex = sapply(seurat_objects, function(s) {
    s$sex
  }) %>% unlist()
) %>%
  arrange(dataset)

p1 <- ggplot(all_qc_data, aes(x = dataset, y = genes, fill = sex)) +
  geom_violin(scale="width", width = 0.6, position = position_dodge(0.7)) +
  xlab("") +
  ylab("Genes / cell") +
  ggtitle("Genes / cell") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

p2 <- ggplot(all_qc_data, aes(x = dataset, y = counts, fill = sex)) +
  geom_violin(scale="width", width = 0.6, position = position_dodge(0.7)) +
  xlab("") +
  ylab("Counts / cell") +
  ggtitle("Counts / cell") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

p3 <- ggplot(all_qc_data, aes(x = dataset, y = mt_perc, fill = sex)) +
  geom_violin(scale="width", width = 0.6, position = position_dodge(0.7)) +
  xlab("") +
  ylab("% mitochondrial genes") +
  ylim(0, 100) +
  ggtitle("% mitochondrial genes") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

if (plot_out_type == "pdf") {
  pdf("plots/all_qc.pdf", width = 15, height = 5)
} else if (plot_out_type == "png") {
  png("plots/all_qc.png", width = 1500, height = 500)
}

grid.arrange(p1, p2, p3, ncol = 3)

if (plot_out_type != "none") {
  dev.off()
}

# Scale with SCT -> PCA -> UMAP
seurat_objects_SCT <- sapply(seurat_objects, function(s) {
  print(paste("Processing", s$author[1], s$year[1], "-", s$sex[1]))

  s %>%
    SCTransform() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    RunPCA(npcs=50) %>%
    RunUMAP(dims = 1:20) -> seuratobj

  outfile <- paste0("rds_outs/", s$orig.ident[1], "_SCT.rds")
  saveRDS(seuratobj, file = outfile)

  return(seuratobj)
})

# Load filed, if needed
# filenames <- list.files("rds_outs/", pattern = "_SCT.rds")
# seurat_objects_SCT <- pblapply(filenames, readRDS)

allQC <- lapply(seurat_objects_SCT, function(s) {
  data.frame(
    dataset = paste(s$author[1], s$year[1], "-", s$sex[1]),
    sex = factor(s$sex[1], levels = c("M", "F")),
    avg_gene_cell = mean(s$nFeature_RNA),
    avg_counts_cell = mean(s$nCount_RNA),
    avg_mt_perc = mean(s$percent_mt)
  )
})

if (plot_out_type == "pdf") {
  pdf("plots/qc_summary_all_datasets.pdf",
    width = 15, height = 10,
    pointsize = 25
  )
} else if (plot_out_type == "png") {
  png("plots/qc_summary_all_datasets.png",
    width = 1500, height = 1000,
    pointsize = 25
  )
}

do.call(rbind, allQC) %>%
  pivot_longer(
    cols = -c(dataset, sex),
    names_to = "qc_feature",
    values_to = "value"
  ) %>%
  # Rename feature names to be more readable
  mutate(qc_feature = case_when(
    qc_feature == "avg_gene_cell" ~ "Avg genes/cell",
    qc_feature == "avg_counts_cell" ~ "Avg counts/cell",
    qc_feature == "avg_mt_perc" ~ "Avg % mitochondrial genes"
  )) %>%
  # Reorder levels
  mutate(qc_feature = factor(qc_feature, levels = c(
    "Avg genes/cell",
    "Avg counts/cell",
    "Avg % mitochondrial genes"
  ))) %>%
  ggplot(aes(y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(x = 0, y = value), size = 2) +
  ylab("") +
  xlab("") +
  # Extend y scale to 0
  scale_y_continuous(limits = c(0, NA)) +
  xlim(-1, 1) +
  facet_wrap(sex ~ qc_feature, scales = "free_y") +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    title = element_text(size = 16),
    strip.text = element_text(size = 15)
  )

if (plot_out_type != "none") {
  dev.off()
}

if (plot_out_type == "pdf") {
  pdf("plots/umap_all_datasets.pdf",
    width = 15, height = 10,
    units = "in", res = 300
  )
} else if (plot_out_type == "png") {
png("plots/umap_all_datasets.png",
  width = 15, height = 10,
  units = "in", res = 300
)

p <- lapply(seurat_objects_SCT, function(obj) {
  DimPlot(obj, pt.size = 0.05, cols = "lightgray") +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(paste(obj$author[1], obj$year[1], "-", obj$sex[1])) +
    NoLegend() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14),
      title = element_text(size = 15)
    )
})
grid.arrange(grobs = p, ncol = 4)

if (plot_out_type != "none") {
  dev.off()
}