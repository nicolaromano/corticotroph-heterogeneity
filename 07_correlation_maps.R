# 07_correlation_maps.R
# Calculates correlation maps between the markers of subclusters of the datasets

library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(grid)
library(gridExtra)
library(scales)
library(pbapply)
library(pheatmap)
library(nlme)
library(emmeans)
library(pbapply)

# Set the pseudo-random seed for reproducibility, changing this will change
# the random swaps, however it should not affect the results in any significant way
set.seed(12345)

plot_out_type <- "pdf" # or "pdf" or "none"

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Data to process ("M" or "F")
data_to_process <- "F"

# Number of random swaps to perform
num_swaps <- 100

filenames <- dir("rds_outs",
  pattern = paste0(data_to_process, "_subclustered.rds"),
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)

# Check if the plots/correlation_maps/ directory exists and create it if not
if (!dir.exists("plots/correlation_maps")) {
  dir.create("plots/correlation_maps")
}

get_marker_correlation <- function(reference_id,
                                   corr_type = "pearson",
                                   mrk_logfc_thr = 0.25, mrk_min_pct = 0.2,
                                   n_rnd = 100,
                                   do_plot = TRUE) {
  #' Create correlation maps between the markers of subclusters of one dataset
  #' and all other dataset
  #' @param reference_id: the id of the reference study
  #' @param corr_type: the method for calculating the correlation (or the cosine similarity);
  #' one of "pearson" (default), "kendall", "spearman", or "cosine"
  #' @param mrk_logfc_thr: the minimum increase in log FC to consider - defaults to 0.25
  #' @param mrk_min_pct: the minimum % of cells expressing the marker - defaults to 0.2
  #' @param n_rnd: number of random swaps to perform
  #' @param do_plot: whether to plots the maps (defaults to TRUE) and the average correlation
  #' @return a list containing the correlation of the markers of each subcluster of the
  #' reference dataset either between themself (on a per-dataset/per-subcluster basis, expected to be high);
  #' between markers of one cluster vs the other clusters (expected to be low); as well as the markers
  #' vs themselves in two randomized versions of the matrix, one where we randomly swap gene expression between
  #' cells (expected to be ~0) and one where we swap gene expression only between cells of the same subcluster

  stopifnot(corr_type %in% c("pearson", "kendall", "spearman", "cosine"))

  dataset_ids <- as.character(sapply(seurat_corticotrophs, function(obj) {
    paste(obj$author[1], obj$year[1], "-", obj$sex[1])
  }))

  # Find marker genes for our reference dataset
  markers <- FindAllMarkers(seurat_corticotrophs[[reference_id]],
    logfc.threshold = mrk_logfc_thr, min.pct = mrk_min_pct,
    only.pos = TRUE
  )

  markers %>%
    # Get only significant markers
    subset(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    # Arrange by p value
    # IMPORTANT - from the manual
    # Unlike other dplyr verbs, arrange() largely ignores grouping;
    # you need to explicitly mention grouping variables (or use .by_group = TRUE)
    # in order to group by them
    arrange(-p_val_adj, .by_group = TRUE) %>%
    select(cluster, gene, p_val_adj) %>%
    # Get the top 50
    top_n(n = 50, wt = -p_val_adj) -> markers_filtered

  corr_matrices <- lapply(seq_along(seurat_corticotrophs), function(i) {
    print(paste("PROCESSING", dataset_ids[i]))
    # Get the expression matrix, and filter it for the marker genes of the reference
    # (or at least the ones that are there!)
    # gene_filter <- intersect(markers_filtered$gene, rownames(seurat_corticotrophs[[i]]))

    print("Calculating correlation")
    expr <- as.matrix(GetAssayData(seurat_corticotrophs[[i]])[markers_filtered$gene, ])
    # Calculate the correlation matrix
    if (corr_type == "cosine") {
      # Cosine similarity
      # From https://stats.stackexchange.com/questions/31565
      expr_cor <- expr / sqrt(rowSums(expr * expr))
      expr_cor <- expr_cor %*% t(expr_cor)
    } else {
      expr_cor <- cor(t(expr), method = corr_type)
    }
    # This is necessary as sometimes the same marker can appear in different clusters
    genenames <- make.unique(colnames(expr_cor))
    colnames(expr_cor) <- genenames
    rownames(expr_cor) <- genenames

    # Randomly swap gene expression between cells (we do NOT swap between genes
    # to keep expression distribution for each gene)
    print("Generating random matrices - global")
    expr_rnd_cor <- pblapply(1:n_rnd, function(x) {
      expr_rnd <- t(apply(expr, 1, sample))
      if (corr_type == "cosine") {
        expr_rnd_cor <- expr_rnd / sqrt(rowSums(expr_rnd * expr_rnd))
        expr_rnd_cor <- expr_rnd_cor %*% t(expr_rnd_cor)
      } else {
        expr_rnd_cor <- cor(t(expr_rnd), method = corr_type)
      }

      colnames(expr_rnd_cor) <- genenames
      rownames(expr_rnd_cor) <- genenames
      return(expr_rnd_cor)
    })

    # Randomly swap gene expression between cells of each cluster
    # (we do NOT swap between genes to keep expression distribution
    # for each gene, nor between clusters to check cluster robustness)
    print("Generating random matrices - per subcluster")
    expr_rnd_within_cor <- pblapply(1:n_rnd, function(x) {
      expr_rnd_within <- expr
      for (cl in unique(Idents(seurat_corticotrophs[[i]])))
      {
        cl_barcodes <- which(Idents(seurat_corticotrophs[[i]]) == cl)
        expr_rnd_within[, cl_barcodes] <-
          t(apply(expr_rnd_within[, cl_barcodes], 1, sample))
      }

      if (corr_type == "cosine") {
        expr_rnd_within_cor <- expr_rnd_within / sqrt(rowSums(expr_rnd_within * expr_rnd_within))
        expr_rnd_within_cor <- expr_rnd_within_cor %*% t(expr_rnd_within_cor)
      } else {
        expr_rnd_within_cor <- cor(t(expr_rnd_within), method = corr_type)
      }

      colnames(expr_rnd_within_cor) <- genenames
      rownames(expr_rnd_within_cor) <- genenames
      return(expr_rnd_within_cor)
    })

    return(list(
      cor = expr_cor,
      rnd_cor = expr_rnd_cor,
      rnd_within_cor = expr_rnd_within_cor
    ))
  })

  names(corr_matrices) <- paste0(dataset_ids, "_vs_", dataset_ids[reference_id])

  if (do_plot) {
    cluster_palette <- hue_pal()(length(unique(markers_filtered$cluster)))
    names(cluster_palette) <- as.character(unique(markers_filtered$cluster))
    corr_plots <- lapply(seq_along(corr_matrices), function(i) {
      m <- corr_matrices[[i]]
      # We might miss some of the gene names in some of the correlation matrices
      # Getting gene names each time here helps avoiding mismatches
      clusters <- data.frame(cluster = markers_filtered$cluster)
      rownames(clusters) <- colnames(m$cor)
      clusters$cluster <- as.character(clusters$cluster)

      # Positions of the gaps in the heatmap, between clusters
      # This automatically checks where the cluster changes
      gaps <- which(diff(as.numeric(clusters$cluster)) == 1)

      dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
      dataset_names <- unlist(dataset_names)
      main_title <- paste(
        "Markers from", dataset_names[2],
        "\ndata from", dataset_names[1]
      )

      # from the pheatmap help page:
      # br: a sequence of numbers that covers the range of values in mat
      # **and is one element longer than color vector**.
      # Used for mapping values to colors.
      br <- seq(-0.5, 0.5, length.out = 257)
      anno_colors <- list(cluster = cluster_palette[seq_along(unique(clusters$cluster))])

      p <- pheatmap(m$cor,
        show_rownames = FALSE, show_colnames = FALSE,
        annotation_row = clusters,
        annotation_col = clusters,
        color = hcl.colors(256, "PRGn"),
        na_col = "white",
        cluster_rows = FALSE, cluster_cols = FALSE,
        gaps_row = gaps, gaps_col = gaps,
        main = main_title,
        breaks = br,
        border_color = NA,
        legend = FALSE,
        annotation_colors = anno_colors,
        fontsize = 10,
        silent = TRUE
      )

      return(p$gtable)
    })

    # Make an heatmat with the same color scale
    ph <- pheatmap(
      cor(
        matrix(rnorm(100), ncol = 10),
        matrix(rnorm(100), ncol = 10)
      ),
      color = hcl.colors(256, "PRGn"),
      breaks = seq(-0.3, 0.3, length.out = 257),
      silent = TRUE
    )
    # Grab the legend
    leg <- ph$gtable$grobs[[4]]

    # Append it to the plot list
    corr_plots_with_leg <- append(corr_plots, list(leg), after = 3)

    png(paste0("plots/correlation_maps/heatmap_markers_", dataset_ids[reference_id], ".png"),
      width = 15, height = 10, units = "in", res = 300
    )

    grid.arrange(
      grobs = corr_plots_with_leg, ncol = 4,
      widths = c(1, 1, 1, 0.3), top = textGrob(paste(
        "Correlation of marker genes - Reference dataset:",
        dataset_ids[reference_id]
      ), gp = gpar(fontsize = 16, fontface = "bold"))
    )
    dev.off()

    # Do the same, but for the randomly shuffled data.
    # We select a random run (from the n_rnd that we ran) to show

    run_id <- sample(1:n_rnd, size = 1)
    corr_plots_rnd <- lapply(seq_along(corr_matrices), function(i) {
      m <- corr_matrices[[i]]
      # We might miss some of the gene names in some of the correlation matrices
      # Getting gene names each time here helps avoiding mismatches
      clusters <- data.frame(cluster = markers_filtered$cluster)
      rownames(clusters) <- colnames(m$rnd_cor[[run_id]])
      clusters$cluster <- as.character(clusters$cluster)

      # Positions of the gaps in the heatmap, between clusters
      # This automatically checks where the cluster changes
      gaps <- which(diff(as.numeric(clusters$cluster)) == 1)

      dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
      dataset_names <- unlist(dataset_names)
      main_title <- paste(
        "Markers from", dataset_names[2], "\ndata from",
        dataset_names[1], "(shuffled)"
      )
      br <- seq(-0.5, 0.5, length.out = 257)
      anno_colors <- list(cluster = cluster_palette[seq_along(unique(clusters$cluster))])

      p <- pheatmap(m$rnd_cor[[run_id]],
        show_rownames = FALSE, show_colnames = FALSE,
        annotation_row = clusters,
        annotation_col = clusters,
        color = hcl.colors(256, "PRGn"),
        na_col = "white",
        cluster_rows = FALSE, cluster_cols = FALSE,
        gaps_row = gaps, gaps_col = gaps,
        main = main_title,
        breaks = br,
        border_color = NA,
        legend = FALSE,
        annotation_colors = anno_colors,
        fontsize = 8,
        silent = TRUE
      )

      return(p$gtable)
    })

    corr_plots_rnd_with_leg <- append(corr_plots_rnd, list(leg), after = 3)

    png(paste0("plots/correlation_maps/heatmap_markers_rnd_", dataset_ids[reference_id], ".png"),
      width = 15, height = 10, units = "in", res = 300
    )
    grid.arrange(
      grobs = corr_plots_rnd_with_leg, ncol = 4,
      widths = c(1, 1, 1, 0.3), top = textGrob(paste(
        "Correlation of marker genes - Reference dataset:",
        dataset_ids[reference_id], "(shuffled)"
      ), gp = gpar(fontsize = 16, fontface = "bold"))
    )
    dev.off()

    corr_plots_rnd_within <- lapply(seq_along(corr_matrices), function(i) {
      m <- corr_matrices[[i]]
      # We might miss some of the gene names in some of the correlation matrices
      # Getting gene names each time here helps avoiding mismatches
      clusters <- data.frame(cluster = markers_filtered$cluster[match(
        colnames(m$rnd_within_cor[[run_id]]),
        markers_filtered$gene
      )])
      rownames(clusters) <- colnames(m$rnd_within_cor[[run_id]])
      clusters$cluster <- as.character(clusters$cluster)

      # Positions of the gaps in the heatmap, between clusters
      # This automatically checks where the cluster changes
      gaps <- which(diff(as.numeric(clusters$cluster)) == 1)

      dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
      dataset_names <- unlist(dataset_names)
      main_title <- paste(
        "Markers from", dataset_names[2], "\ndata from",
        dataset_names[1], "(shuffled in clust)"
      )
      br <- seq(-0.5, 0.5, length.out = 257)
      anno_colors <- list(cluster = cluster_palette[seq_along(unique(clusters$cluster))])

      p <- pheatmap(m$rnd_within_cor[[run_id]],
        show_rownames = FALSE, show_colnames = FALSE,
        annotation_row = clusters,
        annotation_col = clusters,
        color = hcl.colors(256, "PRGn"),
        na_col = "white",
        cluster_rows = FALSE, cluster_cols = FALSE,
        gaps_row = gaps, gaps_col = gaps,
        main = main_title,
        breaks = br,
        border_color = NA,
        legend = FALSE,
        annotation_colors = anno_colors,
        fontsize = 8,
        silent = TRUE
      )

      return(p$gtable)
    })

    corr_plots_rnd_within_with_leg <- append(corr_plots_rnd_within, list(leg), after = 3)

    png(paste0("plots/correlation_maps/heatmap_markers_rnd_within_", dataset_ids[reference_id], ".png"),
      width = 15, height = 10, units = "in", res = 300
    )

    grid.arrange(
      grobs = corr_plots_rnd_within_with_leg,
      ncol = 4, widths = c(1, 1, 1, 0.3), top = ""
    )
    dev.off()
  } # End of if(do_plot)

  return(list(markers = markers_filtered, corr = corr_matrices))
}

get_average_correlations <- function(reference_id, corr_type = "pearson", do_plot = TRUE) {
  #' Gets the average correlation matrix for within/between clusters
  #' and randomly shuffled matrices
  #' @param reference_id: the id of the reference dataset
  #' @param do_plot: whether to do plots, defaults to TRUE
  #' @return a dataframe with the average correlation values
  res <- get_marker_correlation(
    reference_id = reference_id,
    corr_type = corr_type,
    n_rnd = num_swaps,
    do_plot = do_plot
  )

  markers <- res$markers
  corr_matr <- res$corr

  list_names <- names(corr_matr)
  corr_matr <- lapply(seq_along(corr_matr), function(i) {
    m <- corr_matr[[i]]

    # Set diagonals of correlation matrix to NA
    diag(m$cor) <- NA
    m$rnd_cor <- lapply(m$rnd_cor, function(m) {
      diag(m) <- NA
      return(m)
    })
    m$rnd_within_cor <- lapply(m$rnd_within_cor, function(m) {
      diag(m) <- NA
      return(m)
    })

    datasets <- unlist(strsplit(names(corr_matr)[i], split = "_vs_"))
    if (datasets[1] == datasets[2]) {
      m$group <- "within"
    } else {
      m$group <- "between"
    }

    m$mean_cor <- markers %>%
      group_by(cluster) %>%
      group_map(function(genes, cluster) {
        genes <- unlist(genes$gene)
        # Ensure we haven't dropped any gene
        genes <- genes[genes %in% rownames(m$cor)]
        other_genes <- markers$gene[!markers$gene %in% genes]
        other_genes <- other_genes[other_genes %in% rownames(m$cor)]

        list(
          real_within = mean(m$cor[genes, genes], na.rm = TRUE),
          real_between = mean(m$cor[genes, other_genes], na.rm = TRUE),
          rnd = mean(
            sapply(m$rnd_cor, function(mtx) {
              mean(mtx[genes, genes], na.rm = TRUE)
            }) # end sapply
          ), # end mean
          rnd_within = mean(
            sapply(m$rnd_within_cor, function(mtx) {
              mean(mtx[genes, genes], na.rm = TRUE)
            }) # end sapply
          ) # end mean
        ) # end list
      }) # end group_map

    return(m)
  })

  names(corr_matr) <- list_names

  grp_names <- c(
    "Within cluster", "Between clusters",
    "Random shuffle", "Random shuffle within"
  )

  # Calculate average correlation
  mean_cor <- sapply(corr_matr, FUN = function(m) {
    mean_cor <- data.frame(
      Correlation = unlist(m$mean_cor),
      Group = factor(grp_names, levels = grp_names)
    )

    return(mean_cor)
  }, simplify = FALSE, USE.NAMES = TRUE)


  mean_cor %>%
    imap(function(data, dataset) {
      data %>%
        mutate("Dataset" = dataset) %>%
        separate(Dataset, sep = "_vs_", into = c("Dataset1", "Dataset2")) %>%
        mutate("SameDataset" = factor(Dataset1 == Dataset2, levels = c(TRUE, FALSE)))
    }) -> corr_matr_all

  corr_matr_all <- do.call("rbind", c(corr_matr_all, make.row.names = FALSE))

  # We don't care about random shuffling in different datasets
  corr_matr_all <- corr_matr_all %>%
    subset(!(Group %in% c("Random shuffle", "Random shuffle within") & SameDataset == FALSE))

  if (do_plot) {
    if (plot_out_type == "pdf") {
      pdf(paste0("plots/correlation_maps/avg_corr_", data_to_process, "_", reference_id, ".pdf"),
        width = 15, height = 8
      )
    } else if (plot_out_type == "png") {
      png(paste0("plots/correlation_maps/avg_corr_", data_to_process, "_", reference_id, ".png"),
        width = 15, height = 8, units = "in", res = 300
      )
    }

    ggplot(corr_matr_all, aes(Group, Correlation)) +
      geom_boxplot(outlier.shape = NA, aes(fill = SameDataset)) +
      scale_fill_manual(values = c("#86A93F", "#B56492"), name = "Same dataset") +
      geom_point(aes(col = SameDataset), position = position_jitterdodge(jitter.width = 0.1)) +
      scale_color_manual(values = c("#253700", "#6F1749"), name = "Same dataset") +
      scale_x_discrete(labels = label_wrap(10)) +
      ylim(c(-0.1, 0.3)) +
      theme(
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12)
      )

    if (plot_out_type != "none") {
      dev.off()
    }
  }

  return(corr_matr_all)
}

all_avg_corr <- sapply(seq_along(seurat_corticotrophs),
  FUN = get_average_correlations, do_plot = FALSE, corr_type = "pearson",
  simplify = FALSE,
  USE.NAMES = TRUE
)

all_avg_corr <- do.call("rbind", all_avg_corr)

write.csv(all_avg_corr, file = paste0("all_average_correlations", data_to_process, ".csv"), row.names = FALSE)

all_avg_corr <- read.csv(paste0("all_average_correlations", data_to_process, ".csv"))

if (plot_out_type == "pdf") {
  pdf(paste0("plots/correlation_maps/all_avg_corr_", data_to_process, ".pdf"),
    width = 15, height = 8
  )
} else if (plot_out_type == "png") {
  png(paste0("plots/correlation_maps/all_avg_corr_", data_to_process, ".png"),
    width = 15, height = 8, units = "in", res = 300
  )
}

all_avg_corr %>%
  subset(Group %in% c("Within cluster", "Between clusters")) %>%
  ggplot(aes(Group, Correlation)) +
  geom_boxplot(outlier.shape = NA, aes(fill = SameDataset)) +
  scale_fill_manual(values = c("#86A93F", "#B56492"), name = "Same dataset") +
  geom_point(aes(col = SameDataset), position = position_jitterdodge(jitter.width = 0.1)) +
  scale_color_manual(values = c("#253700", "#6F1749"), name = "Same dataset") +
  scale_x_discrete(labels = label_wrap(10)) +
  geom_hline(yintercept = 0, lty = "dotted") +
  ylim(c(-0.1, 0.3)) +
  theme(
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap(~Dataset1) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

if (plot_out_type != "none") {
  dev.off()
}

plot_correlation_boxplot <- function(data, with_shuffles = FALSE, ncol = 4) {
  #' Plots a boxplot of the correlation values
  #' @param data: the data to plot  
  #' @param with_shuffles: whether to include the random shuffles in the plot
  #' @return NULL, but saves the plot
  
  if (!with_shuffles) {
    data <- data %>%
      subset(Group %in% c("Within cluster", "Between clusters"))
  }

  ggplot(data, aes(Group, Correlation)) +
    geom_boxplot(
      outlier.shape = NA, aes(fill = SameDataset),
      position = position_dodge2(preserve = "single")
    ) +
    scale_fill_manual(values = c("#86A93F", "#B56492"), name = "Same dataset") +
    geom_point(aes(col = SameDataset),
      position = position_jitterdodge(jitter.width = 0.1)
    ) +
    scale_color_manual(values = c("#253700", "#6F1749"), name = "Same dataset") +
    scale_x_discrete(labels = label_wrap(10)) +
    geom_hline(yintercept = 0, lty = "dotted") +
    ylim(c(-0.1, 0.3)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14)
    ) +
    facet_wrap(~Dataset1, ncol = ncol)
}

if (plot_out_type == "pdf") {
  pdf(paste0("plots/correlation_maps/all_avg_corr", data_to_process, ".pdf"),
    width = 15, height = 8
  )
} else if (plot_out_type == "png") {
  png(paste0("plots/correlation_maps/all_avg_corr", data_to_process, ".png"),
    width = 15, height = 8, units = "in", res = 300
  )
}

plot_correlation_boxplot(all_avg_corr, TRUE, ncol=3)

if (plot_out_type != "none") {
  dev.off()
}

# Put M and F together
all_avg_corr_mf <- rbind(
  read.csv("all_average_correlationsM.csv"),
  read.csv("all_average_correlationsF.csv")
)

if (plot_out_type == "pdf") {
    pdf("plots/correlation_maps/all_avg_corr_MF.pdf", width = 9, height = 8)
  } else if (plot_out_type == "png") {
    png("plots/correlation_maps/all_avg_corr_MF.png", width = 9, height = 8, units = "in", res = 300)
}

plot_correlation_boxplot(all_avg_corr_mf, with_shuffles = FALSE, ncol=3)

if (plot_out_type != "none") {
  dev.off()
}

model <- all_avg_corr_mf %>%
  mutate(Sex = str_sub(Dataset1, -1)) %>%
  drop_na() %>%
  lme(Correlation ~ Group + SameDataset + Sex,
    random = ~ 1 | Dataset1, data = .
  )

summary(model)
intervals(model)[[1]] %>% 
  as.data.frame() %>% 
  format(digits = 2, scientific = FALSE)
# pairwise comparisons
emmeans(model, pairwise ~ Group | SameDataset) %>% confint
emmeans(model, pairwise ~ SameDataset | Group) %>% confint

# Get mean corr of random shuffle
all_avg_corr_mf %>% 
  dplyr::filter(Group == "Random shuffle") %>% 
  dplyr::group_by(Dataset1) %>%
  dplyr::summarise(mean_cor = mean(Correlation, na.rm = TRUE)) -> mean_rnd

mean(mean_rnd$mean_cor)
sd(mean_rnd$mean_cor)
