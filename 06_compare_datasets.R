# 06_compare_datasets.R
# Compare the datasets based on their common markers

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(tibble)
library(pbapply)
library(pheatmap)
library(igraph)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

####### Compare datasets #######

# Data to process ("M" or "F")
data_to_process <- "F"
output_format <- "pdf" # "pdf", "png" or "none"

markers_output_folder <- "markers/"
if (!dir.exists(markers_output_folder)) {
  dir.create(markers_output_folder)
  print(paste("Created output directory", markers_output_folder))
}

filenames <- dir("rds_outs",
  pattern = paste0(data_to_process, "_subclustered.rds"),
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)
names(seurat_corticotrophs) <- sapply(seurat_corticotrophs, function(s) {
  s$orig.ident[1]
})

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
for (i in seq_along(all_markers)) {
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
  png(paste0("plots/common_markers_boxplot_", data_to_process, ".png"),
    width = 5, height = 10,
    units = "in", res = 300
  )
} else if (output_format == "pdf") {
  pdf(paste0("plots/common_markers_boxplot_", data_to_process, ".pdf"),
    width = 10, height = 5
  )
}

# For each cluster in each dataset, get the maximum % of common markers
common_markers %>%
  # Marker correspondance is transitive, but we want
  # to plot the maximum percentage both ways
  # i.e. if A and B have 10% common markers then
  # we want to plot 10% for A-B and 10% for B-A
  bind_rows(
    common_markers %>%
      mutate(
        Tmp = Dataset2_name,
        Dataset2_name = Dataset1_name,
        Dataset1_name = Tmp,
        Tmp = Dataset2,
        Dataset2 = Dataset1,
        Dataset1 = Tmp,
        Tmp = Cluster2,
        Cluster2 = Cluster1,
        Cluster1 = Tmp
      )
  ) %>%
  select(-Tmp) %>%
  group_by(Dataset1_name, Dataset2_name, Cluster1) %>%
  summarise(Max_Pct = max(Percentage)) %>%
  ungroup() %>%
  mutate(Dataset = str_split(Dataset1_name, " ", simplify = TRUE)[, 1]) %>%
  mutate(Dataset2 = str_split(Dataset2_name, " ", simplify = TRUE)[, 1]) %>%
  ggplot(aes(x = Dataset, y = Max_Pct)) +
  geom_boxplot(aes(fill = Dataset2),
    outlier.size = 0.7,
    outlier.alpha = 1, alpha = 0.9
  ) +
  scale_fill_brewer(
    palette = ifelse(data_to_process == "M", "Greens",
      "Purples"), name = "Dataset"
  ) +
  # Global box
  geom_boxplot(
    # Note these works because All becomes the first factor
    # as is the first in alphabetical order...
    # I could not find a better way to do this
    data = . %>% mutate(Dataset = factor("All")),
    aes(x = Dataset, y = Max_Pct), fill = "lightgray"
  ) +
  # geom_point(aes(col= Dataset2), size = 2,
  #   position = position_dodge(width = 0.75),
  # ) +
  ylim(0, 100) +
  ylab("Maximum percentage of common markers") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  )

if (output_format != "none") {
  dev.off()
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
  width <- 14
} else {
  width <- 8
}

if (output_format == "png") {
  png(paste0("plots/perc_common_markers_", data_to_process, ".png"),
    width = width, height = 7, units = "in", res = 300
  )
} else if (output_format == "pdf") {
  pdf(paste0("plots/perc_common_markers_", data_to_process, ".pdf"),
    width = width, height = 7
  )
}

grid.arrange(grobs = pl_list, ncol = 4)

if (output_format != "none") {
  dev.off()
}

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
      l <- in_radius * layout_in_circle(subgraph(network, nodes_by_group[[i]]$FullName))
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
  } else if (output_format == "png") {
    png(paste0("plots/common_markers_graph_thr", thr, "_", data_to_process, ".png"),
      width = 7, height = 7, units = "in", res = 300
    )
  }

  if (data_to_process == "M") {
    community_palette <- c(
      "#F16745", "#FFC65D",
      "#7BC8A4", "#C44A9A",
      "#4A90F2", "#8A5A44"
    )
  } else {
    community_palette <- c(
      "#F11E4E", "#5B2C6F",
      "#1A9337", "#4CC3D9"
    )
  }
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
  memberships <- membership(cluster_walktrap(graph, steps = 100))

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
    # "Lone" communities of only 1 node should be marked as NA
    seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]][tb[seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]]] == 1] <- NA

    # Convert to factor
    seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]] <- factor(seurat_corticotrophs[[i]]@meta.data[[paste0("marker_community_", thr)]], levels = sort(names(tb[tb > 1])))
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
  } else if (output_format == "png") {
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

  if (output_format != "none") {
    dev.off()
  }

  community_prop <- lapply(seurat_corticotrophs, function(s) {
    table(Community = s$marker_community_10) %>%
      as.data.frame() %>%
      filter(Freq > 0) %>%
      mutate(Freq = Freq / ncol(s) * 100) %>%
      rbind(data.frame(Community = "None", Freq = 100 - sum(.$Freq))) %>%
      mutate(Dataset = paste(s$author[1]))
  })


  if (output_format == "pdf") {
    pdf(paste0("plots/markers_community_proportions_thr", thr, "_", data_to_process, ".pdf"),
      width = 7, height = 7
    )
  } else if (output_format == "png") {
    png(paste0("plots/markers_community_proportions_", data_to_process, ".png"),
      width = 7, height = 7, units = "in", res = 300
    )
  }

  community_prop <- do.call(rbind, community_prop)

  p <- ggplot(community_prop, aes(x = Dataset, y = Freq, fill = fct_rev(Community))) +
    scale_fill_manual(values = c(comm_colors, c("None" = "lightgray")), name = "Community") +
    geom_bar(stat = "identity") +
    ylab("Percentage of cells") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  print(p)

  if (output_format != "none") {
    dev.off()
  }
}

# filenames <- dir("rds_outs",
#   pattern = paste0(data_to_process, "_subclustered.rds"),
#   full.names = TRUE
# )
# seurat_corticotrophs <- pblapply(filenames, readRDS)

# Community-specific markers
for (perc_similarity in c(10, 15, 20)) {
  community_mark_intersection <- list()
  community_mark_union <- list()
  num_clusters_per_community <- list()

  for (i in seq_along(seurat_corticotrophs)) {
    obj <- seurat_corticotrophs[[i]]
    dataset_name <- names(seurat_corticotrophs)[i]
    for (cl in unique(Idents(obj))) {
      comm_num <- obj[[paste0("marker_community_", perc_similarity)]][Idents(obj) == cl, ][1]

      if (is.na(comm_num)) {
        next
      }

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
  comm_names <- comm_names[!is.na(comm_names)]

  # Write to file
  for (cn in comm_names) {
    community_un <- community_mark_union[[cn]]
    community_in <- community_mark_intersection[[cn]]

    print(paste(
      "Community", cn, " - markers intersection / union = ",
      format(length(community_in) / length(community_un) * 100, digits = 2), "%"
    ))

    write.table(community_un, file.path(markers_output_folder, paste0(cn, "_", data_to_process, "_perc_simil_", perc_similarity, "_markers_union.csv")),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    write.table(community_in, file.path(markers_output_folder, paste0(cn, "_", data_to_process, "_perc_simil_", perc_similarity, "_markers_intersection.csv")),
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
        markers_output_folder, paste0("markers_communities_", perc_similarity, "_", data_to_process, ".csv")
      ),
      row.names = FALSE
    )
} # end for perc_similarity

marker_stats <- data.frame()

for (simil in c(10, 15, 20)) {
  markers_filenames <- dir(markers_output_folder,
    pattern = paste0(data_to_process, "_perc_simil_", simil, ".*_union.csv"),
    full.names = TRUE
  )
  total_markers <- sapply(markers_filenames, function(f) {
    read.table(f) %>%
      nrow()
  })

  markers_filenames <- dir(markers_output_folder,
    pattern = paste0(data_to_process, "_perc_simil_", simil, ".*_intersection.csv"),
    full.names = TRUE
  )
  common_markers <- sapply(markers_filenames, function(f) {
    # Need to use try as some files are empty and read.table throws an error
    marker_data <- NULL
    try(marker_data <- read.table(f), silent = TRUE)
    if (is.null(marker_data)) {
      return(0)
    }

    return(nrow(marker_data))
  })

  communities <- sapply(strsplit(markers_filenames, "_"), function(x) {
    x[2]
  })
  marker_stats <- rbind(
    marker_stats,
    data.frame(
      Similarity = simil,
      Community = communities,
      Total_Markers = total_markers,
      Common_Markers = common_markers
    )
  )
} # end for simil

write.csv(marker_stats, file.path(markers_output_folder, paste0("markers_stats_", data_to_process, ".csv")), row.names = FALSE)

community_df <- read.csv(file.path(markers_output_folder, paste0("markers_communities_", 10, "_", data_to_process, ".csv")))

plot_perc_cell_per_marker <- function(perc_similarity = 10) {
  #' Plots the percentage of cells in each dataset, for each community that express common markers
  #' @param perc_similarity: optional, defaults to 10. The similarity threshold we're using

  files <- list.files(markers_output_folder, pattern = paste0(data_to_process, "_perc_simil_", perc_similarity, ".*_union.csv"), full.names = TRUE)
  community_summary <- read.csv(paste0(markers_output_folder, "markers_communities_", perc_similarity, "_", data_to_process, ".csv")) %>%
    subset(!is.na(Community))

  markers <- lapply(files, function(f) {
    df <- data.frame(read.table(f)) %>%
      rename(Gene = 1) %>%
      mutate(Community = paste("Community", strsplit(f, "_")[[1]][[2]]))
  })

  markers <- do.call(rbind, markers)

  marker_expression <- lapply(seurat_corticotrophs, function(obj) {
    expr_per_cluster <- sapply(unique(Idents(obj)), function(id) {
      # Get expression for each cluster

      expr <- GetAssayData(subset(obj, ident = id))[markers$Gene, ]
      # mean the expression only in cells that express the gene...
      expr <- apply(expr, 1, function(x) {
        x <- x[x > 0]
        (x - mean(x)) / sd(x)
      })
      expr <- sapply(expr, mean)
      # perc <- apply(expr, 1, function(x){length(x[x>0]) / length(x) * 100})
    })

    colnames(expr_per_cluster) <- factor(unique(Idents(obj)))

    expr_long <- expr_per_cluster %>%
      as.data.frame() %>%
      rownames_to_column("Gene") %>%
      pivot_longer(cols = -Gene, values_to = "Expression", names_to = "Cluster") %>%
      mutate(Cluster = as.integer(sub("V", "", Cluster))) %>%
      left_join(markers, by = "Gene", relationship = "many-to-many") %>%
      subset(!is.na(Community))

    return(expr_long)
  })

  marker_expression <- lapply(names(marker_expression), function(x) {
    marker_expression[[x]] <- marker_expression[[x]] %>%
      mutate(Dataset = x)
  })

  marker_expression <- do.call(rbind, marker_expression)

  marker_expression %>%
    subset(!is.na(Expression)) %>%
    group_by(Dataset, Cluster, Community) %>%
    summarize(Mean_Expression = mean(Expression)) %>%
    ggplot(aes(x = Community, y = Cluster)) +
    geom_tile(aes(fill = Mean_Expression)) +
    scale_fill_gradient2(low = "#00aeff", mid = "#ffffe4", high = "#e96363") +
    facet_grid(Dataset ~ 1, scales = "free") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 5))
}

# Save the seurat objects back to disk
pblapply(seurat_corticotrophs, function(s) {
  saveRDS(s, paste0("rds_outs/", s$orig.ident[1], "_subclustered.rds"))
})

####### GENERAL DATASET STATS #######

# Load all, M+F
filenames <- dir("rds_outs",
  pattern = paste0("_subclustered.rds"),
  full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)
names(seurat_corticotrophs) <- sapply(seurat_corticotrophs, function(s) {
  s$orig.ident[1]
})

# How much of the reads are Pomc?
pomc_perc <- lapply(seurat_corticotrophs, function(s) {
  expr <- GetAssayData(s, slot = "counts")
  sum(expr["Pomc", ]) / sum(expr) * 100
})

quantile(unlist(pomc_perc))

# How many reads are 0?
zero_perc <- lapply(seurat_corticotrophs, function(s) {
  expr <- GetAssayData(s, slot = "counts")
  sum(expr == 0) / (nrow(expr) * ncol(expr)) * 100
})

quantile(unlist(zero_perc))