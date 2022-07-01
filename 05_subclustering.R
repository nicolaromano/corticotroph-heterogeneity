library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pbapply)
library(pheatmap)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Read only male files
filenames <- dir("rds_outs/", pattern = "M_corticotrophs.rds")
filenames <- paste0("rds_outs/", filenames)

# Ho 2020 has only a handful of cells - we'are not considering it
filenames <- filenames[-grep("Ho2020", filenames)]
seurat_corticotrophs <- pblapply(filenames, readRDS)

##### CLUSTERING QUALITY ASSESSMENT #####

### Elbow method ###
# We look at the "elbow" in a plot of resolution vs WSS
# (Within-Cluster-Sum of Squared Errors)
elbowMethod <- function(obj, resolutions = seq(0, 2, 0.1))
  {
  #' Calculated the WSS at various clustering resolutions
  #' This performs clustering on the data at various resolutions
  #' @param obj: the Seurat object
  #' @param resol: the resolutions to test (default is seq(0, 2, 0.1))
  WSS <- NULL
  for (resol in resolutions)
    {
    # The results of this clustering get stored in:
    # 1) obj$seurat_clusters (which is a shortcut for obj@metadata$seurat_clusters)
    # 2) obj$SCT_snn_res.xxx where xxx is the resolution
    # 1) is more convenient to use in this function, but mind that it gets
    # overwritten every time FindClusters is run!
    obj <- FindClusters(obj, resolution = resol)
    
    # Get UMAP coordinates
    umap_coords <- data.frame(obj@reductions$umap@cell.embeddings)
    
    # Find the centroids of each cluster
    umap_coords %>% 
      group_by(cluster = obj$seurat_clusters) %>% 
      summarise(UMAP_1 = mean(UMAP_1),
                UMAP_2 = mean(UMAP_2)) %>% 
      ungroup() -> centroids
    centroids$cluster <- as.integer(centroids$cluster)
    # Find the Within-Cluster-Sum of Squared Errors (WSS) 
    # This is the sum of the squared distance of each point 
    # in a cluster to its centroid
    umap_coords <- umap_coords %>% 
      mutate(cluster = as.integer(obj$seurat_clusters))
    
    res <- sum(apply(umap_coords, 1, function(row)
    {
      dist(rbind(row[1:2], centroids[centroids$cluster == row[3], 2:3])) ^ 2
    }))
    
    WSS <- rbind(WSS, c(resol = resol, 
                        WSS = res, 
                        n_clusters = length(levels(obj$seurat_clusters))))
  }
  
  WSS <- data.frame(WSS)
  
  ggplot(WSS, aes(x = resol, y = WSS)) +
    geom_line() +
    geom_point() +
    geom_text(aes(label = n_clusters), vjust = -1) +
    xlab("Resolution") +
    # Works best with uniformly spaced (and not too many) resolutions
    scale_x_continuous(breaks = seq(min(resolutions), max(resolutions), 
                                    length.out = length(resolutions))) + 
    ggtitle(obj$orig.ident[1]) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
  }


elbow_plots <- lapply(seurat_corticotrophs, function(obj)
  {
  obj <- obj %>% 
    RunUMAP(dims = 1:20) %>% 
    FindNeighbors

  elbowMethod(obj, resolutions = seq(0.1, 1.5, 0.2))
  })

do.call("grid.arrange", c(elbow_plots, ncol=3))

# These were manually chosen
resolutions <- c(0.3, 0.3, 0.5, 0.5, 0.3, 0.3)

seurat_corticotrophs <- lapply(seq_along(seurat_corticotrophs), function(i) {
  obj <- seurat_corticotrophs[[i]]
  obj <- obj %>% 
    RunUMAP(dims = 1:20) %>% 
    FindNeighbors %>% 
    FindClusters(resolution = resolutions[i])
  
  return(obj)
  })

cort_plots <- lapply(seurat_corticotrophs, function(obj){
  DimPlot(obj) +
    xlab(expression(UMAP[1])) +
    ylab(expression(UMAP[2])) +
    ggtitle(obj$orig.ident[1])
})  

png("plots/corticotrophs_all_datasets.png", width = 800, height = 1200)
do.call("grid.arrange", cort_plots)
dev.off()

# Save to RDS
pblapply(seurat_corticotrophs, function(obj){
  saveRDS(obj, file = paste0("rds_outs/", obj$orig.ident[1], "_subclustered.rds"))
  })

##### Find marker genes #####
i <- 4
obj <- seurat_corticotrophs[[i]]
markers <- FindAllMarkers(obj, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE)
markers %>% 
  subset(p_val_adj < 0.05) %>% 
  group_by(cluster) %>%
  arrange(-p_val_adj) %>% 
  select(gene) %>% 
  top_n(50) -> markers_filtered
  
expr <- as.matrix(GetAssayData(obj)[markers$gene,])
expr_cor <- cor(t(expr))
expr_rnd_cor <- cor(matrix(sample(expr), ncol = nrow(expr), nrow = ncol(expr)))
genenames <- make.unique(colnames(expr_cor))
colnames(expr_cor) <- genenames
rownames(expr_cor) <- genenames

colnames(expr_rnd_cor) <- genenames
rownames(expr_rnd_cor) <- genenames

clusters <- data.frame(cluster = markers$cluster)
rownames(clusters) <- genenames

pheatmap(expr_cor, show_rownames = FALSE, show_colnames = FALSE, 
         annotation_row = clusters, 
         na_col = "white",
         cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(expr_rnd_cor, show_rownames = FALSE, show_colnames = FALSE, 
         annotation_row = clusters, 
         # annotation_col = clusters,
         cluster_rows = FALSE, cluster_cols = FALSE)

mean_cor <- markers_filtered %>% 
  group_by(cluster) %>% 
  group_map(function(genes, cluster){
    genes <- unlist(genes)
    other_genes <- markers_filtered$gene[!markers_filtered$gene %in% genes]
    list(real = mean(expr_cor[genes,genes], na.rm = TRUE),
         rnd = mean(expr_rnd_cor[genes,genes], na.rm = TRUE))
    })

mean_cor <- matrix(unlist(mean_cor), ncol = 2)
colnames(mean_cor) <- c("Real", "Random")

boxplot(mean_cor)
