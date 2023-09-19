library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(pbapply)
library(igraph)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

# Data to process ("M" or "F")
data_to_process <- "M"
output_format <- "pdf" # or "png"

markers_output_folder <- "markers"

filenames <- dir("rds_outs",
  pattern = paste0(data_to_process, "_subclustered.rds"),
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)
names(seurat_corticotrophs) <- gsub("(\\.rds|rds_outs/|_subclustered)", "", filenames)

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

  markers %>%
    filter(p_val_adj < 0.05)
})

# Save all markers to csv

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

for (i in 1:length(all_markers)) {
  write.csv(all_markers[[i]], file.path(markers_output_folder, paste0(names(all_markers)[i], "_markers.csv")))
}

# Find common markers
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

if (output_format == "png") {
  png("plots/percentage_common_markers.png",
    width = 5, height = 10,
    units = "in", res = 300
  )
} else {
  pdf("plots/percentage_common_markers.pdf",
    width = 5, height = 10
  )
}

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

if (output_format == "png") {
  png(paste0("plots/perc_common_markers_", data_to_process, ".png"),
    width = width, height = 7, units = "in", res = 300
  )
} else {
  pdf(paste0("plots/perc_common_markers_", data_to_process, ".pdf"),
    width = width, height = 7
  )
}

grid.arrange(grobs = pl_list, ncol = 4)
dev.off()
##### Similarity graph ######

format_dataset_name <- function(name) {
  #' Go from Xxxxxxxxx-n to Xn
  #' @param name: The name to format
  #' @return: The formatted name
  #' @examples
  #' format_dataset_name("Lopez2021M-2") will return 'L2'
  #' format_dataset_name("Vennekens2021M-10") will return 'V10'
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

  if (output_format == "pdf") {
    pdf(paste0("plots/common_markers_graph_thr", thr, "_", data_to_process, ".pdf"),
      width = 7, height = 7
    )
  } else {
    png(paste0("plots/common_markers_graph_thr", thr, "_", data_to_process, ".png"),
      width = 7, height = 7, units = "in", res = 300
    )
  }

  community_palette <- c(
    "#F16745", "#FFC65D", "#7BC8A4",
    "#4CC3D9", "#93648D", "#808070"
  )

  get_similarity_graph(thr,
    out_radius = out_radius, in_radius = in_radius,
    node_palette = community_palette
  )
  dev.off()
}

for (thr in c(10, 15, 20)) {
  graph <- get_similarity_graph(thr,
    do_plot = FALSE
  )
  memberships <- membership(walktrap.community(graph, steps = 100))

  tb <- table(memberships)

  # Assign a colour to each community.
  comm_colors <- community_palette[sort(unique(memberships))]
  # Communities with a single node are grey
  comm_colors[tb == 1] <- "lightgrey"
  names(comm_colors) <- names(tb)

  for (i in seq_along(seurat_corticotrophs)) {
    # Add a column to the metadata with the community
    seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]] <-
      memberships[match(
        paste0(
          substr(seurat_corticotrophs[[i]]$author[1], 1, 1),
          Idents(seurat_corticotrophs[[i]])
        ),
        names(memberships)
      )]
    seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]] <- factor(seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]], levels = sort(unique(memberships)))
  }

  # Now plot the UMAP reductions, colouring by community
  community_plots <- lapply(seurat_corticotrophs, function(s) {
    # Get the first letter of the dataset name.
    # This corresponds to the names of the graph nodes
    graph_nodes_name <- substr(s$author[1], 1, 1)

    palette <- comm_colors

    p <- DimPlot(s, group.by = paste0("marker_community_", thr)) +
      scale_color_manual(
        name = "Community",
        values = palette
      ) +
      ggtitle(paste(s$author[1], s$year[1], "-", s$sex[1])) +
      theme(legend.position = "none")
  })

  if (output_format == "pdf") {
    pdf(paste0("plots/umap_communities_common_markers_thr", thr, "_", data_to_process, ".pdf"),
      width = 7, height = 7
    )
  } else {
    png(paste0("plots/umap_communities_common_markers_thr", thr, "_", data_to_process, ".png"),
      width = 7, height = 7, units = "in", res = 300
    )
  }

  grid.arrange(
    grobs = community_plots, ncol = 3,
    top = textGrob(
      paste("Cells with >= ", thr, "% common markers"),
      gp = gpar(fontsize = 20)
    )
  )
  dev.off()
}

# Save the seurat objects back to disk
pblapply(seurat_corticotrophs, function(s) {
  saveRDS(s, paste0("rds_outs/", s$orig.ident[1], "_subclustered.rds"))
})

# Community-specific markers
community_mark_intersection <- list()
community_mark_union <- list()
num_clusters_per_community <- list()

perc_similarity <- 15

for (i in seq_along(seurat_corticotrophs)) {
  obj <- seurat_corticotrophs[[i]]
  dataset_name <- names(seurat_corticotrophs)[i]
  for (cl in unique(Idents(obj))) {
    comm_num <- obj[[paste0("marker_community_", perc_similarity)]][Idents(obj) == cl, ][1]

    # INTERSECTION
    if (is.null(community_mark_intersection[[paste0("comm_", comm_num)]])) {
      community_mark_intersection[[paste0("comm_", comm_num)]] <- all_markers[[dataset_name]] %>%
        filter(cluster == cl) %>%
        pull(gene)
    } else {
      community_mark_intersection[[paste0("comm_", comm_num)]] <- intersect(
        community_mark_intersection[[paste0("comm_", comm_num)]],
        all_markers[[dataset_name]] %>%
          filter(cluster == cl) %>%
          pull(gene)
      )
    }

    # UNION
    if (is.null(community_mark_union[[paste0("comm_", comm_num)]])) {
      community_mark_union[[paste0("comm_", comm_num)]] <- all_markers[[dataset_name]] %>%
        filter(cluster == cl) %>%
        pull(gene)
    } else {
      community_mark_union[[paste0("comm_", comm_num)]] <- union(
        community_mark_union[[paste0("comm_", comm_num)]],
        all_markers[[dataset_name]] %>%
          filter(cluster == cl) %>%
          pull(gene)
      )
    }

    # Number of clusters per community
    if (is.null(num_clusters_per_community[[paste0("comm_", comm_num)]])) {
      num_clusters_per_community[[paste0("comm_", comm_num)]] <- 1
    } else {
      num_clusters_per_community[[paste0("comm_", comm_num)]] <- num_clusters_per_community[[paste0("comm_", comm_num)]] + 1
    }
  } # end for cl
} # end for i

num_clusters_per_community <- unlist(num_clusters_per_community)
# Only get the communities with more than one cluster
comm_names <- sort(names(num_clusters_per_community[num_clusters_per_community > 1]))

# Write to file
for (cn in comm_names) {
  community_un <- community_mark_union[[cn]]
  community_in <- community_mark_intersection[[cn]]

  print(paste(
    "Community", cn, " - markers intersection / union = ",
    format(length(community_in) / length(community_un) * 100, digits = 2), "%"
  ))

  write.table(community_un, file.path(markers_output_folder, paste0(cn, "_markers_union.csv")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  write.table(community_in, file.path(markers_output_folder, paste0(cn, "_markers_intersection.csv")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

# Write a table of Dataset, Cluster, Community
community_df <- lapply(seq_along(seurat_corticotrophs), function(i) {
  obj <- seurat_corticotrophs[[i]]
  dataset_name <- names(seurat_corticotrophs)[i]

  clusters <- unique(Idents(obj))
  sapply(clusters, function(cl) {
    subset(obj, idents = cl)[[paste0(
      "marker_community_",
      perc_similarity
    )]][1, 1]
  }) -> communities

  res <- data.frame(
    Dataset = dataset_name,
    Cluster = clusters,
    Community = communities
  )
})

do.call("rbind", community_df) %>%
  write.csv(
    file.path(
      markers_output_folder, paste0("markers_communities_", perc_similarity, ".csv")
    ),
    row.names = FALSE
  )

plot_gene_rankings <- function(seurat_corticotrophs, reference_dataset, outfile = NULL) {
  #' Plots the gene rankings for each dataset, using one
  #' dataset as a reference.
  #'
  #' Params
  #' ------
  #' seurat_corticotrophs: A list of Seurat objects.
  #' reference_dataset: The dataset to use as a reference (as an index in the list)
  #' outfile: The file to save the plot to. If NULL, the plot is not saved.

  # Get the gene rankings for each dataset
  rankings <- lapply(seurat_corticotrophs, function(s) {
    expr <- GetAssayData(s, slot = "data")
    data.frame(
      Dataset = paste(s$author[1], s$year[1], "-", s$sex[1]),
      Gene = rownames(s),
      Rank = rank(-rowMeans(expr)),
      AvgExpr = rowMeans(expr)
    ) %>%
      arrange(Rank)
  })

  # Get the ranked gene name for the reference dataset
  ref_rankings <- rankings[[reference_dataset]] %>%
    arrange(Rank) %>%
    pull(Gene)

  # Now sort the rankings by the reference dataset
  ranks <- lapply(rankings, function(r) {
    r[ref_rankings, ]$Rank
  })

  if (!is.null(outfile)) {
    png(outfile, width = 15, height = 10, units = "in", res = 300)
  }

  g <- do.call("rbind", rankings) %>%
    as.data.frame() %>%
    # Substitute the ranks with the ranks sorted by the reference dataset's order
    mutate(Rank = Reduce(c, ranks)) %>%
    ggplot(aes(x = Rank, AvgExpr)) +
    geom_line() +
    scale_y_log10() +
    ylab("Average expression") +
    facet_wrap(~Dataset) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 16)
    )

  print(g)

  if (!is.null(outfile)) {
    dev.off()
  }
}

plot_gene_rankings(seurat_corticotrophs, 1, paste0("plots/gene_rankings_", data_to_process, ".png"))
