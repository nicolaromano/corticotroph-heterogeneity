library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(pbapply)
library(igraph)

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
  markers <- FindAllMarkers(d,
    logfc.threshold = 0.25,
    min.pct = 0.2, only.pos = TRUE
  )
})

common_markers <- data.frame(
  Dataset1_name = character(),
  Dataset2_name = character(),
  Dataset1 = numeric(), Dataset2 = numeric(),
  Cluster1 = numeric(), Cluster2 = numeric(),
  Percentage = numeric()
)

for (d1 in 1:(length(seurat_corticotrophs) - 1)) {
  for (d2 in (d1 + 1):length(seurat_corticotrophs)) {
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

ggplot(common_markers, aes(y = Percentage)) +
  geom_boxplot()

prettify_df_name <- function(name) {
  name <- as.character(name)
  name <- substr(name, 1, nchar(name) - 1) # Remove sex
  paste(
    substr(name, 1, nchar(name) - 4), # Author
    substr(name, nchar(name) - 3, nchar(name))
  ) # Year
}

# Get all of the possible combinations of datasets
pl_list <- apply(unique(common_markers[, c("Dataset1", "Dataset2")]), 1, function(d1d2) {
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
    theme(
      axis.title = element_text(size = 10),
      legend.position = "none"
    )
  g
})

g <- ggplot(data.frame(x = 8, y = 0), aes(x, y)) +
  geom_point(aes(size = x, col = x > 10)) +
  scale_size(
    breaks = c(10, 20, 30), limits = c(0, 30), name = element_blank(),
    guide = guide_legend(override.aes = list(colour = rgb(0.5, 0.6, 0.8)))
  ) +
  scale_color_manual(
    values = rgb(0, 0, 0, 0.2), labels = "<10%", name = "% common markers",
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  theme_minimal()


leg <- get_legend(g)
leg$layout[2, ]$b <- 2
pl_list$legend <- leg
do.call("grid.arrange", c(pl_list, ncol = 4))

##### Similarity graph ######


format_dataset_name <- function(name) {
  # Go from Xxxxx20xx-n to Xn
  paste0(substr(name, 1, 1), substr(name, nchar(name), nchar(name)))
}

plot_similarity_graph <- function(min_pct, node_palette = NA,
                                  edge_palette = NA,
                                  label_color = "black") {
  #' Plots a graph, connecting subcluster with at least `min_pct`
  #' markers in common
  #' @param min_pct: The minimum percentage of markers in common
  #' @param node_palette: The palette for the nodes - NA for default (Set2)
  #' @param edge_palette: The palette for the edges - NA for default (Accent)
  #' @param label_color: The colour for the labels - default is black

  common_markers %>%
    mutate(
      From = paste(Dataset1_name, Cluster1, sep = "-"),
      To = paste(Dataset2_name, Cluster2, sep = "-")
    ) %>%
    mutate(
      From = format_dataset_name(From),
      To = format_dataset_name(To)
    ) %>%
    select(c(From, To, Percentage)) %>%
    filter(Percentage >= min_pct) -> markers_nodes

  # We find all of the combinations of dataset name and subcluster number
  datasets <- unique(c(
    paste(common_markers$Dataset1_name,
      common_markers$Cluster1,
      sep = "-"
    ),
    paste(common_markers$Dataset2_name,
      common_markers$Cluster2,
      sep = "-"
    )
  ))
  # Simplify the naming (for graphical purposes) and sort alphabetically
  datasets <- sort(format_dataset_name(datasets), decreasing = TRUE)

  # Now create a graph using the subcluster as vertices and connecting those
  # with >= 10% shared markers
  network <- graph_from_data_frame(
    d = markers_nodes, vertices = datasets,
    directed = FALSE
  )

  # Create a layout for our network

  data.frame(FullName = datasets) %>%
    mutate(
      Dataset = substr(FullName, 1, 1),
      SubCluster = substr(FullName, 2, 2)
    ) %>%
    group_split(Dataset) -> nodes_by_group

  # Find communities
  wc <- walktrap.community(network)
  node_comm <- membership(wc)

  if (is.na(node_palette)) {
    node_palette <- brewer.pal(length(unique(substr(datasets, 1, 1))), "Set2")
  }

  if (is.na(edge_palette)) {
    edge_palette <- brewer.pal(sum(sapply(communities(wc), length) > 1), "Accent")
  }

  # We put the nodes from the same datasets in a circle, then put
  # each dataset on a larger circle
  out_radius <- 5
  in_radius <- 1

  layout <- sapply(seq_along(nodes_by_group), function(i) {
    l <- in_radius * layout.circle(subgraph(network, nodes_by_group[[i]]$FullName))
    l[, 1] <- l[, 1] + out_radius * sin(2 * pi * i / length(nodes_by_group))
    l[, 2] <- l[, 2] + out_radius * cos(2 * pi * i / length(nodes_by_group))

    l
  })

  # Combine all of the sub-layouts into a big layout!
  layout <- do.call("rbind", rev(layout))

  # Define colors for the nodes and edges, based on the dataset of origin
  node_colors <- node_palette[as.integer(factor(substr(V(network)$name, 1, 1)))]
  edge_colors <- edge_palette[membership(wc)[get.edgelist(network)[, 1]]]

  # Finally plot the graph
  plot(network,
    vertex.color = node_colors,
    edge.width = E(network)$Percentage / 5, # Edge width
    edge.curved = 0.3,
    edge.color = edge_colors,
    vertex.label.family = "Arial",
    vertex.label.color = label_color,
    layout = layout
  )
  # layout = layout.fruchterman.reingold)
}


plot_similarity_graph(10, node_palette <- c(
  "#F16745", "#FFC65D", "#7BC8A4",
  "#4CC3D9", "#93648D", "#808070"
),
label_color = rep(c("black", "white", "black"), c(12, 9, 9))
)
