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

# Data to process ("M" or "F")
data_to_process <- "M"

# Read only male files
filenames <- dir("rds_outs",
  pattern = paste0(data_to_process, "_subclustered.rds"),
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)
names(seurat_corticotrophs) <- gsub("(\\.rds|rds_outs/)", "", filenames)

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
          Dataset1_name = paste(
            seurat_corticotrophs[[d1]]$author[1],
            seurat_corticotrophs[[d1]]$year[1], "-",
            seurat_corticotrophs[[d1]]$sex[1]
          ),
          Dataset2_name = paste(
            seurat_corticotrophs[[d2]]$author[1],
            seurat_corticotrophs[[d2]]$year[1], "-",
            seurat_corticotrophs[[d2]]$sex[1]
          ),
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

png("plots/percentage_common_markers.png",
  width = 5, height = 10,
  units = "in", res = 300
)
ggplot(common_markers, aes(y = Percentage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(x = 0), width = 0.1, size = 2) +
  xlim(-0.7, 0.7) +
  ylab("Percentage of common markers") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  )
dev.off()

# Get all of the possible combinations of datasets
pl_list <- apply(unique(common_markers[, c("Dataset1", "Dataset2")]), 1, function(d1d2) {
  common_markers %>%
    subset(Dataset1 == d1d2["Dataset1"]) %>%
    subset(Dataset2 == d1d2["Dataset2"]) -> sub_markers


  g <- ggplot(sub_markers, aes(x = Cluster1, y = Cluster2)) +
    geom_point(aes(size = Percentage, col = Percentage > 10)) +
    scale_size(breaks = c(10, 20, 30), limits = c(0, 30)) +
    scale_color_manual(values = c(rgb(0, 0, 0, 0.2), rgb(0.6, 0.6, 0.8))) +
    xlab(sub_markers$Dataset1_name[1]) +
    ylab(sub_markers$Dataset2_name[1]) +
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

if (data_to_process == "F") {
  width <- 9
} else {
  width <- 8
}

png(paste0("plots/perc_common_markers_", data_to_process, ".png"),
  width = width, height = 7, units = "in", res = 300
)
grid.arrange(grobs = pl_list, ncol = 4)
dev.off()
##### Similarity graph ######

format_dataset_name <- function(name) {
  # Go from Xxxxx20xx-n to Xn
  paste0(substr(name, 1, 1), substr(name, nchar(name), nchar(name)))
}

get_similarity_graph <- function(min_pct, do_plot = TRUE,
                                 node_palette = NULL,
                                 label_color = "black",
                                 out_radius = 5, in_radius = 1) {
  #' Plots a graph, connecting subcluster with at least `min_pct`
  #' markers in common
  #' @param min_pct: The minimum percentage of markers in common
  #' @param do_plot: Whether to plot the graph - default is TRUE
  #' @param node_palette: The palette for the nodes - NULL for default (Set2)
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
  # with >= min_pct% shared markers
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
  wc <- walktrap.community(network, steps = 100)
  node_comm <- membership(wc)
  nodes_per_comm <- table(node_comm)

  if (do_plot) {
    if (!length(node_palette)) {
      node_palette <- brewer.pal(length(unique(substr(datasets, 1, 1))), "Set2")
    }

    node_colors <- node_palette[node_comm]
    # Communities with a single node are white
    node_colors[nodes_per_comm[node_comm] == 1] <- "white"

    # We put the nodes from the same datasets in a circle, then put
    # each dataset on a larger circle

    layout <- sapply(seq_along(nodes_by_group), function(i) {
      l <- in_radius * layout.circle(subgraph(network, nodes_by_group[[i]]$FullName))
      l[, 1] <- l[, 1] + out_radius * sin(2 * pi * i / length(nodes_by_group))
      l[, 2] <- l[, 2] + out_radius * cos(2 * pi * i / length(nodes_by_group))

      l
    })

    # Combine all of the sub-layouts into a big layout!
    layout <- do.call("rbind", rev(layout))

    # Finally plot the graph
    plot(network,
      vertex.color = node_colors,
      edge.width = E(network)$Percentage / 5, # Edge width
      edge.curved = 0.3,
      edge.color = rgb(0, 0, 0, 0.2),
      vertex.label.family = "Arial",
      vertex.label.color = label_color,
      layout = layout,
      # layout = layout.fruchterman.reingold
      main = paste0("Clusters with >= ", min_pct, "% common markers")
    )
  }

  # Return the graph object
  return(network)
}

for (thr in c(10, 15, 20)) {
  if (data_to_process == "M") {
    out_radius <- 5
    in_radius <- 1
  } else {
    out_radius <- 5
    in_radius <- 1.5
  }

  png(paste0("plots/graph_common_clusters_thr", thr, "_", data_to_process, ".png"),
    width = 7, height = 7, units = "in", res = 300
  )
  get_similarity_graph(thr,
    out_radius = out_radius, in_radius = in_radius,
    node_palette = c(
      "#F16745", "#FFC65D", "#7BC8A4",
      "#4CC3D9", "#93648D", "#808070"
    )
  )
  dev.off()
}

for (thr in c(10, 15, 20)) {
  graph <- get_similarity_graph(thr,
    do_plot = FALSE
  )
  memberships <- membership(walktrap.community(graph, steps = 100))

  tb <- table(memberships)

  community_palette <- c(
    "#F16745", "#FFC65D", "#7BC8A4",
    "#4CC3D9", "#93648D", "#808070"
  )

  comm_colors <- community_palette[sort(unique(memberships))]
  comm_colors[tb == 1] <- "lightgrey"
  names(comm_colors) <- names(tb)

  for (i in seq_along(seurat_corticotrophs)) {
    # Add a column to the metadata with the community
    seurat_corticotrophs[[i]][[paste0("marker_community_", thr)]] <-
      memberships[match(paste0(
        substr(seurat_corticotrophs[[i]]$author[1], 1, 1),
        Idents(seurat_corticotrophs[[i]])
      ), names(memberships))]
    seurat_corticotrophs[[i]][[paste0("marker_community_", thr)]] <- factor(seurat_corticotrophs[[i]][[paste0("marker_community_", thr)]], levels = sort(unique(memberships)))
  }

  # Now plot the UMAP reductions, colouring by community
  community_plots <- lapply(seurat_corticotrophs, function(s) {
    # Get the first letter of the dataset name.
    # This corresponds to the names of the graph nodes
    graph_nodes_name <- substr(s$author[1], 1, 1)

    palette <- comm_colors[unique(s$community)]

    p <- DimPlot(s, group.by = paste0("marker_community_", thr)) +
      scale_color_manual(
        name = "Community",
        values = palette
      ) +
      ggtitle(paste(s$author[1], s$year[1], "-", s$sex[1])) +
      theme(legend.position = "none")
  })

  png(paste0("plots/umap_common_clusters_thr", thr, "_", data_to_process, ".png"),
    width = 7, height = 7, units = "in", res = 300
  )
  grid.arrange(
    grobs = community_plots, ncol = 3,
    top = textGrob(
      paste("Cells with >= ", thr, "% common markers"),
      gp = gpar(fontsize = 20)
    )
  )
  dev.off()
}