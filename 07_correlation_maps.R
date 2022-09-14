library(Seurat)
library(ggplot2)
library(gridExtra)
library(scales)
library(tidyverse)
library(pbapply)
library(pheatmap)
library(nlme)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Read only male files
filenames <- dir("rds_outs", pattern = "M_subclustered.rds", 
                 full.names = TRUE)

seurat_corticotrophs <- pblapply(filenames, readRDS)

prettify_df_name <- function(name)
  {
  name = as.character(name)
  name <- substr(name, 1, nchar(name)-1) # Remove sex
  paste(substr(name, 1, nchar(name)-4), # Author
        substr(name, nchar(name)-3, nchar(name))) # Year
  }

get_marker_correlation <- function(reference_id,
                                    corr_type = "pearson",
                                    mrk_logfc_thr = 0.25, mrk_min_pct = 0.2,
                                    n_rnd = 100,
                                    do_plot = TRUE)
  {
  #' Create correlation maps between the markers of subclusters of one dataset
  #' and all other dataset
  #' @param reference_id: the id of the reference study
  #' @param corr_type: the method for calculating the correlation; 
  #' one of "pearson" (default), "kendall", or "spearman"
  #' @param mrk_logfc_thr: the minimum increase in log FC to consider - defaults to 0.25
  #' @param mrk_min_pct: the minimum % of cells expressing the marker - defaults to 0.2
  #' @param n_rnd: number of random swaps to perform
  #' @param do_plot: whether to plots the maps (defaults to TRUE) and the average correlation
  #' @return a list containing the correlation of the markers of each subcluster of the 
  #' reference dataset either between themself (on a per-dataset/per-subcluster basis, expected to be high); 
  #' between markers of one cluster vs the other clusters (expected to be low); as well as the markers
  #' vs themselves in two randomized versions of the matrix, one where we randomly swap gene expression between
  #' cells (expected to be ~0) and one where we swap gene expression only between cells of the same subcluster
  
  stopifnot(corr_type %in% c("pearson", "kendall", "spearman"))
  
  dataset_ids <- as.character(sapply(seurat_corticotrophs, function(obj){obj$orig.ident[1]}))
  
  # Find marker genes for our reference dataset
  markers <- FindAllMarkers(seurat_corticotrophs[[reference_id]], 
                            logfc.threshold = mrk_logfc_thr, min.pct = mrk_min_pct, 
                            only.pos = TRUE)
    
  markers %>% 
    # Get only significan markers
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
    top_n(n = 50, wt = p_val_adj) -> markers_filtered
  
    corr_matrices <- lapply(seq_along(seurat_corticotrophs), function(i){
      print(paste("PROCESSING", dataset_ids[i]))
      # Get the expression matrix, and filter it for the marker genes of the reference
      # (or at least the ones that are there!)
      gene_filter = intersect(markers_filtered$gene, rownames(seurat_corticotrophs[[i]]))
      
      print("Calculating correlation")
      expr <- as.matrix(GetAssayData(seurat_corticotrophs[[i]])[gene_filter,])
      # Calculate the correlation matrix
      expr_cor <- cor(t(expr), method = corr_type)
      # This is necessary as sometimes the same marker can appear in different clusters
      genenames <- make.unique(colnames(expr_cor))
      colnames(expr_cor) <- genenames
      rownames(expr_cor) <- genenames

      # Randomly swap gene expression between cells (we do NOT swap between genes
      # to keep expression distribution for each gene)
      print("Generating random matrices - global")
      expr_rnd_cor <- pblapply(1:n_rnd, function(x) {
        expr_rnd <- t(apply(expr, 1, sample))
        expr_rnd_cor <- cor(t(expr_rnd), method = corr_type)
        
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
        
        expr_rnd_within_cor <- cor(t(expr_rnd_within), method = corr_type)
        
        colnames(expr_rnd_within_cor) <- genenames
        rownames(expr_rnd_within_cor) <- genenames
        return(expr_rnd_within_cor)
      })
      
      return(list(cor = expr_cor, 
                  rnd_cor = expr_rnd_cor, 
                  rnd_within_cor = expr_rnd_within_cor))
     })
    
    names(corr_matrices) <- paste0(dataset_ids, "_vs_", dataset_ids[reference_id])
    
    if (do_plot)
      {
      cluster_palette <-c("0"="#f1c7dd", "1"="#966eac", "2"="#e3337e",
                           "3"="#827775", "4"="#7bcbc0", "5"="#f05129",
                           "6"="#b09977", "7"="#b7cc94", "8"="#f5bd42",
                           "9"="#0b326b")
      corr_plots <- lapply(seq_along(corr_matrices), function(i)
        {
        m = corr_matrices[[i]]
        # We might miss some of the gene names in some of the correlation matrices
        # Getting gene names each time here helps avoiding mismatches
        clusters <- data.frame(cluster = markers_filtered$cluster[match(colnames(m$cor), 
                                                                        markers_filtered$gene)])
        rownames(clusters) <- colnames(m$cor)
        clusters$cluster = as.character(clusters$cluster)
        
        # Positions of the gaps in the heatmap, between clusters
        # This automatically checks where the cluster changes
        gaps <- which(diff(as.numeric(clusters$cluster)) == 1)
        
        dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
        dataset_names <- unlist(dataset_names)
        main_title <- paste("Markers from", prettify_df_name(dataset_names[2]), 
                            "\ndata from", prettify_df_name(dataset_names[1]))

        # from the pheatmap help page:
        # br: a sequence of numbers that covers the range of values in mat 
        # **and is one element longer than color vector**. 
        # Used for mapping values to colors. 
        br <- seq(-0.4, 0.4, length.out=257)
        anno_colors <- list(cluster = cluster_palette[1:length(unique(clusters$cluster))])
        
        p <- pheatmap(m$cor, show_rownames = FALSE, show_colnames = FALSE, 
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
                 silent = TRUE)
        
        return(p$gtable)
        })

      # Make an heatmat with the same color scale
      ph <- pheatmap(cor(matrix(rnorm(100), ncol=10), 
                   matrix(rnorm(100), ncol=10)), 
               color = hcl.colors(256, "PRGn"),
               breaks = seq(-0.3, 0.3, length.out=257),
               silent = TRUE)
      # Grab the legend
      leg <- ph$gtable$grobs[[4]]

      # Append it to the plot list
      corr_plots_with_leg <- append(corr_plots, list(leg), after = 3)
      do.call("grid.arrange", c(corr_plots_with_leg, 
                                list(ncol=4, widths = c(1, 1, 1, 0.3)), top=""))
      
      # Do the same, but for the randomly shuffled data.
      # We select a random run (from the n_rnd that we ran) to show

      run_id = sample(1:n_rnd, size = 1)
      corr_plots_rnd <- lapply(seq_along(corr_matrices), function(i)
      {
        m = corr_matrices[[i]]
        # We might miss some of the gene names in some of the correlation matrices
        # Getting gene names each time here helps avoiding mismatches
        clusters <- data.frame(cluster = markers_filtered$cluster[match(colnames(m$rnd_cor[[run_id]]), 
                                                                        markers_filtered$gene)])
        rownames(clusters) <- colnames(m$rnd_cor[[run_id]])
        clusters$cluster = as.character(clusters$cluster)
        
        # Positions of the gaps in the heatmap, between clusters
        # This automatically checks where the cluster changes
        gaps <- which(diff(as.numeric(clusters$cluster)) == 1)
        
        dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
        dataset_names <- unlist(dataset_names)
        main_title <- paste("Markers from", dataset_names[2], "\ndata from", 
                            dataset_names[1], "(shuffled)")
        br <- seq(-0.5, 1, 0.01)
        anno_colors <- list(cluster = cluster_palette[1:length(unique(clusters$cluster))])
        
        p <- pheatmap(m$rnd_cor[[run_id]], show_rownames = FALSE, show_colnames = FALSE, 
                      annotation_row = clusters,
                      annotation_col = clusters, 
                      na_col = "white",
                      cluster_rows = FALSE, cluster_cols = FALSE, 
                      gaps_row = gaps, gaps_col = gaps,
                      main = main_title,
                      breaks = br,
                      border_color = NA,
                      legend = FALSE,
                      annotation_colors = anno_colors,
                      fontsize = 8,
                      silent = TRUE)
        
        return(p$gtable)
      })

      corr_plots_rnd_with_leg <- append(corr_plots_rnd, list(leg), after = 3)
      do.call("grid.arrange", c(corr_plots_rnd_with_leg, 
                                list(ncol=4, widths = c(1, 1, 1, 0.3)), top=""))
      
      corr_plots_rnd_within <- lapply(seq_along(corr_matrices), function(i)
      {
        m = corr_matrices[[i]]
        # We might miss some of the gene names in some of the correlation matrices
        # Getting gene names each time here helps avoiding mismatches
        clusters <- data.frame(cluster = markers_filtered$cluster[match(colnames(m$rnd_within_cor[[run_id]]), 
                                                                        markers_filtered$gene)])
        rownames(clusters) <- colnames(m$rnd_within_cor[[run_id]])
        clusters$cluster = as.character(clusters$cluster)
        
        # Positions of the gaps in the heatmap, between clusters
        # This automatically checks where the cluster changes
        gaps <- which(diff(as.numeric(clusters$cluster)) == 1)
        
        dataset_names <- strsplit(names(corr_matrices[i]), split = "_vs_")
        dataset_names <- unlist(dataset_names)
        main_title <- paste("Markers from", dataset_names[2], "\ndata from", 
                            dataset_names[1], "(shuffled in clust)")
        br <- seq(-0.5, 1, 0.01)
        anno_colors <- list(cluster = cluster_palette[1:length(unique(clusters$cluster))])
        
        p <- pheatmap(m$rnd_within_cor[[run_id]], show_rownames = FALSE, show_colnames = FALSE, 
                      annotation_row = clusters,
                      annotation_col = clusters, 
                      na_col = "white",
                      cluster_rows = FALSE, cluster_cols = FALSE, 
                      gaps_row = gaps, gaps_col = gaps,
                      main = main_title,
                      breaks = br,
                      border_color = NA,
                      legend = FALSE,
                      annotation_colors = anno_colors,
                      fontsize = 8,
                      silent = TRUE)
        
        return(p$gtable)
      })
      
      corr_plots_rnd_within_with_leg <- append(corr_plots_rnd_within, list(leg), after = 3)
      do.call("grid.arrange", c(corr_plots_rnd_within_with_leg, 
                                list(ncol=4, widths = c(1, 1, 1, 0.3)), top=""))
    } # if

    return (list(markers = markers_filtered, corr = corr_matrices))
  }

get_average_correlations <- function(reference_id, do_plot = TRUE)
  {
  #' Gets the average correlation matrix for within/between clusters
  #' and randomly shuffled matrices
  #' @param reference_id: the id of the reference dataset
  #' @param do_plot: whether to do plots, defaults to TRUE
  #' @return a dataframe with the average correlation values
  res <- get_marker_correlation(reference_id = reference_id, 
                                do_plot = do_plot)

  markers <- res$markers
  corr_matr <- res$corr
  
  list_names <- names(corr_matr)
  corr_matr <- lapply(seq_along(corr_matr), function(i)
    {
    m <- corr_matr[[i]]
  
    # Set diagonals of correlation matrix to NA  
    diag(m$cor) <- NA
    m$rnd_cor <- lapply(m$rnd_cor, function(m){
      diag(m) <- NA
      return(m)
    })
    m$rnd_within_cor <- lapply(m$rnd_within_cor, function(m){
      diag(m) <- NA
      return(m)
    })
  
    datasets <- unlist(strsplit(names(corr_matr)[i], split = "_vs_"))
    if (datasets[1] == datasets[2])
      m$group = "within"
    else
      m$group = "between"
    
    m$mean_cor <- markers %>%
       group_by(cluster) %>% 
       group_map(function(genes, cluster){
           genes <- unlist(genes$gene)
           # Ensure we haven't dropped any gene
           genes <- genes[genes %in% rownames(m$cor)]
           other_genes <- markers$gene[!markers$gene %in% genes]
           other_genes <- other_genes[other_genes %in% rownames(m$cor)]
           
           list(real_within = mean(m$cor[genes,genes], na.rm = TRUE),
                real_between = mean(m$cor[genes,other_genes], na.rm = TRUE),
                rnd = mean(sapply(m$rnd_cor, function(mtx){
                  mean(mtx[genes,genes], na.rm = TRUE)
                }) # end sapply
                ), # end mean
                rnd_within = mean(sapply(m$rnd_within_cor, function(mtx){
                  mean(mtx[genes,genes], na.rm = TRUE)
                }) # end sapply
                ) # end mean
           ) # end list
         }) # end group_map
    
    return(m)
    })
  
  names(corr_matr) <- list_names
  
  grp_names <- c("Within cluster", "Between clusters",
                 "Random shuffle", "Random shuffle within")
  
  # Calculate average correlation
  mean_cor <- sapply(corr_matr, function(m)
    {
    mean_cor <- data.frame(Correlation = unlist(m$mean_cor),
                             Group = factor(grp_names, levels = grp_names))
    
    return(mean_cor)
    }, simplify = FALSE, USE.NAMES = TRUE)
  
  
  mean_cor %>% 
    imap(function(data, dataset){
      data %>% 
        mutate("Dataset" = dataset) %>% 
        separate(Dataset, sep = "_vs_", into = c("Dataset1", "Dataset2")) %>% 
        mutate("SameDataset" = factor(Dataset1 == Dataset2, levels = c(TRUE, FALSE)))
    }) -> corr_matr_all
    
  corr_matr_all <- do.call("rbind", c(corr_matr_all, make.row.names = FALSE))
  
  # We don't care about random shuffling in different datasets
  corr_matr_all <- corr_matr_all %>% 
    subset(!(Group %in% c("Random shuffle", "Random shuffle within") & SameDataset == FALSE))

  if (do_plot)
    {
    ggplot(corr_matr_all, aes(Group, Correlation)) +
      geom_boxplot(outlier.shape = NA, aes(fill = SameDataset)) +
      scale_fill_manual(values = c("#86A93F", "#B56492"), name = "Same dataset") +
      geom_point(aes(col = SameDataset), position = position_jitterdodge(jitter.width = 0.1)) +
      scale_color_manual(values = c("#253700", "#6F1749"), name = "Same dataset") +
      scale_x_discrete(labels = label_wrap(10)) +
      ylim(c(-0.1, 0.3)) +
      theme(axis.text = element_text(size = 11),
            strip.text = element_text(size = 12))
    }
  
  return(corr_matr_all)
}

all_avg_corr <- sapply(1:length(seurat_corticotrophs), 
                       FUN = get_average_correlations, do_plot=FALSE, 
                       simplify = FALSE,
                       USE.NAMES = TRUE)

all_avg_corr <- do.call("rbind", all_avg_corr)

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
    theme(axis.text = element_text(size = 11),
          strip.text = element_text(size = 12)) +
    facet_wrap(~Dataset1)

ggplot(all_avg_corr, aes(Group, Correlation)) +
  geom_boxplot(outlier.shape = NA, aes(fill = SameDataset), 
               position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#86A93F", "#B56492"), name = "Same dataset") +
  geom_point(aes(col = SameDataset), 
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_color_manual(values = c("#253700", "#6F1749"), name = "Same dataset") +
  scale_x_discrete(labels = label_wrap(10)) +
  geom_hline(yintercept = 0, lty = "dotted") +
  ylim(c(-0.1, 0.3)) +
  theme(axis.text = element_text(size = 11),
        strip.text = element_text(size = 12)) +
  facet_wrap(~Dataset1)

write.csv(all_avg_corr, file = "all_average_correlations.csv", row.names = FALSE)
all_avg_corr %>% 
  drop_na() -> all_avg_corr_nona

model <- lme(Correlation ~ Group + SameDataset, random = ~ 1 | Dataset1, data = all_avg_corr_nona)
summary(model)
