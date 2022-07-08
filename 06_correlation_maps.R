library(Seurat)
library(ggplot2)
library(tidyverse)
library(pbapply)
library(pheatmap)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Read only male files
filenames <- dir("rds_outs", pattern = "M_subclustered.rds", 
                 full.names = TRUE)

seurat_corticotrophs <- pblapply(filenames, readRDS)

##### Find marker genes #####
mean_cor <- lapply(seurat_corticotrophs, function(obj)
  {
  markers <- FindAllMarkers(obj, logfc.threshold = 0.25, min.pct = 0.1, 
                            only.pos = TRUE)
  markers %>% 
    subset(p_val_adj < 0.05) %>% 
    group_by(cluster) %>%
    # IMPORTANT - from the manual
    # Unlike other dplyr verbs, arrange() largely ignores grouping; 
    # you need to explicitly mention grouping variables (or use .by_group = TRUE) 
    # in order to group by them
    arrange(-p_val_adj, .by_group = TRUE) %>% 
    select(cluster, gene, p_val_adj) %>% 
    top_n(n = 50, wt = p_val_adj) -> markers_filtered
  
  
  expr <- as.matrix(GetAssayData(obj)[markers_filtered$gene,])
  expr_cor <- cor(t(expr), method = "pearson")
  
  genenames <- make.unique(colnames(expr_cor))
  colnames(expr_cor) <- genenames
  rownames(expr_cor) <- genenames
  
  n.rnd <- 100
  
  # Randomly swap gene expression between cells (we do NOT swap between genes
  # to keep expression distribution for each gene)
  expr_rnd_cor <- pblapply(1:n.rnd, function(x) {
    
    expr_rnd <- t(apply(expr, 1, sample))
    expr_rnd_cor <- cor(t(expr_rnd), method = "pearson")
    
    colnames(expr_rnd_cor) <- genenames
    rownames(expr_rnd_cor) <- genenames
    return(expr_rnd_cor)
  })
  
  # Randomly swap gene expression between cells of each cluster 
  # (we do NOT swap between genes to keep expression distribution 
  # for each gene, nor between clusters to check cluster robustness)
  expr_rnd_within_cor <- pblapply(1:n.rnd, function(x) {
    expr_rnd_within <- expr
    for (cl in unique(Idents(obj)))
    {
      cl_barcodes <- which(Idents(obj) == cl)
      expr_rnd_within[, cl_barcodes] <- 
        t(apply(expr_rnd_within[, cl_barcodes], 1, sample))
    }
    
    expr_rnd_within_cor <- cor(t(expr_rnd_within), method = "pearson")
    
    colnames(expr_rnd_within_cor) <- genenames
    rownames(expr_rnd_within_cor) <- genenames
    return(expr_rnd_within_cor)
  })
  
  clusters <- data.frame(cluster = markers_filtered$cluster)
  rownames(clusters) <- genenames
  
  # Positions of the gaps in the heatmap, between clusters
  # This automatically checks where the cluster changes
  gaps <- which(diff(as.numeric(clusters$cluster)) == 1)
  
  pheatmap(expr_cor, show_rownames = FALSE, show_colnames = FALSE, 
           annotation_row = clusters,
           annotation_col = clusters, 
           na_col = "white",
           cluster_rows = FALSE, cluster_cols = FALSE, 
           gaps_row = gaps, gaps_col = gaps)
  
  # Plot one of the random correlation matrices
  n <- sample(1:n.rnd, 1)
  pheatmap(expr_rnd_cor[[n]], show_rownames = FALSE, show_colnames = FALSE, 
           na_col = "white",
           annotation_row = clusters, 
           annotation_col = clusters,
           gaps_row = gaps, gaps_col = gaps,
           cluster_rows = FALSE, cluster_cols = FALSE)
  
  pheatmap(expr_rnd_within_cor[[n]], show_rownames = FALSE, show_colnames = FALSE, 
           na_col = "white",
           annotation_row = clusters, 
           annotation_col = clusters,
           gaps_row = gaps, gaps_col = gaps,
           cluster_rows = FALSE, cluster_cols = FALSE)
  
  diag(expr_cor) <- NA
  expr_rnd_cor <- lapply(expr_rnd_cor, function(m){
    diag(m) <- NA
    return(m)
  })
  expr_rnd_within_cor <- lapply(expr_rnd_within_cor, function(m){
    diag(m) <- NA
    return(m)
  })
  
  mean_cor <- markers_filtered %>%
    group_by(cluster) %>% 
    group_map(function(genes, cluster){
      genes <- unlist(genes$gene)
      other_genes <- markers_filtered$gene[!markers_filtered$gene %in% genes]
      
      list(real_within = mean(expr_cor[genes,genes], na.rm = TRUE),
           real_between = mean(expr_cor[genes,other_genes], na.rm = TRUE),
           rnd = mean(sapply(expr_rnd_cor, function(m){
             mean(m[genes,genes], na.rm = TRUE)
           }) # end sapply
           ), # end mean
           rnd_within = mean(sapply(expr_rnd_within_cor, function(m){
             mean(m[genes,genes], na.rm = TRUE)
           }) # end sapply
           ) # end mean
      ) # end list
    }) # end group_map
  
  
  grp_names <- c("Within cluster", "Between clusters", 
                 "Random shuffle", "Random shuffle within")
  mean_cor <- data.frame(Correlation = unlist(mean_cor),
                         Group = factor(grp_names, levels = grp_names))
  
  return(mean_cor)
  })


mean_cor_all <- lapply(seq_along(mean_cor), function(i){
    cbind(mean_cor[[i]], Dataset = i)
  })

mean_cor_all <- do.call("rbind", mean_cor_all)
mean_cor_all$Dataset <- factor(mean_cor_all$Dataset)
ggplot(mean_cor_all, aes(Group, Correlation)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) 
