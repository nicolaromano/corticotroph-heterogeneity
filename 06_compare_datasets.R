library(Seurat)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(dplyr)
library(pbapply)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Read only male files
filenames <- dir("rds_outs",
  pattern = "M_subclustered.rds",
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)

get_common_markers <- function(markers1, markers2, cluster1, cluster2) {
  #' Gets the percentage of marker genes in common between two
  #' clusters from two datasets
  #' @param markers1 - the markers from dataset 1, found using FindAllMarkers
  #' @param markers2 - the markers from dataset 2, found using FindAllMarkers
  #' @param cluster1 - the cluster from dataset 1 to consider
  #' @param cluster2 - the cluster from dataset 2 to consider
  markers1 <- subset(markers1, cluster == cluster1)
  markers2 <- subset(markers2, cluster == cluster2)
  common <- length(intersect(markers1$gene, markers2$gene)) # n genes in common
  total <- length(union(markers1$gene, markers2$gene)) # n genes markers 1
  return(common / total * 100)
}

all_markers <- lapply(seurat_corticotrophs, function(d) {
  markers <- FindAllMarkers(d, logfc.threshold = 0.25, 
                            min.pct = 0.2, only.pos = TRUE)
})

common_markers <- data.frame(
  Dataset1_name = character(),
  Dataset2_name = character(),
  Dataset1 = numeric(), Dataset2 = numeric(),
  Cluster1 = numeric(), Cluster2 = numeric(),
  Percentage = numeric()
)

for (d1 in 1:(length(seurat_corticotrophs)-1)) {
  for (d2 in (d1+1):length(seurat_corticotrophs)) {
    for (cl1 in unique(Idents(seurat_corticotrophs[[d1]]))) {
      for (cl2 in unique(Idents(seurat_corticotrophs[[d2]]))) {
        perc <- get_common_markers(
          all_markers[[d1]], all_markers[[d2]],
          cl1, cl2
        )
        res <- data.frame(
          Dataset1_name = seurat_corticotrophs[[d1]]$orig.ident[1],
          Dataset2_name = seurat_corticotrophs[[d2]]$orig.ident[1],
          Dataset1 = d1, Dataset2 = d2,
          Cluster1 = cl1, Cluster2 = cl2,
          Percentage = perc
        )
        common_markers <- rbind(common_markers, res)
      }
    }
  }
}

# Just to make the data frame prettier
rownames(common_markers) <- NULL

ggplot(common_markers, aes(y=Percentage)) +
  geom_boxplot()

prettify_df_name <- function(name)
  {
  name = as.character(name)
  name <- substr(name, 1, nchar(name)-1) # Remove sex
  paste(substr(name, 1, nchar(name)-4), # Author
        substr(name, nchar(name)-3, nchar(name))) # Year
  }

# Get all of the possible combinations of datasets
pl_list <- apply(unique(common_markers[,c('Dataset1', 'Dataset2')]), 1, function(d1d2) {
  common_markers %>% 
    subset(Dataset1 == d1d2["Dataset1"]) %>% 
    subset(Dataset2 == d1d2["Dataset2"]) -> sub_markers
  
  
  g <- ggplot(sub_markers, aes(x = Cluster1, y = Cluster2)) +
    geom_point(aes(size = Percentage, col = Percentage > 10)) +
    scale_size(breaks = c(10, 20, 30), limits = c(0, 30)) +
    scale_color_manual(values = c(rgb(0, 0, 0, 0.2), rgb(0.6, 0.6, 0.8))) +
    xlab(prettify_df_name(sub_markers$Dataset1_name[1])) +
    ylab(prettify_df_name(sub_markers$Dataset2_name[1])) +
    theme_minimal() +
    theme(axis.title = element_text(size=10),
          legend.position = "none")
  g
  })

g <- ggplot(data.frame(x=8, y=0), aes(x, y)) +
  geom_point(aes(size=x, col = x>10)) +
  scale_size(breaks = c(10, 20, 30), limits = c(0, 30), name = element_blank(), 
             guide = guide_legend(override.aes = list(colour = rgb(0.5, 0.6, 0.8)))) +
  scale_color_manual(values=rgb(0, 0, 0, 0.2), labels = "<10%", name="% common markers",
                     guide = guide_legend(override.aes = list(size = 3))) +
  theme_minimal()


leg <- get_legend(g)
leg$layout[2,]$b = 2
pl_list$legend <- leg
do.call("grid.arrange", c(pl_list, ncol=4))
