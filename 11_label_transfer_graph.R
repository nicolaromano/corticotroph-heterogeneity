library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(pbapply)
library(dplyr)
library(igraph)
library(RColorBrewer)

# Load corticotrophs
datasets <- read.csv("datasets.csv")

data_to_process <- "M" # M or F

datasets <- datasets %>%
    dplyr::filter(author != "Ho") %>%
    dplyr::filter(sex == data_to_process)

output_format <- "pdf" # or "png"
output_folder <- "label_transfer_model_output/predictions/"

filenames <- dir("rds_outs",
    pattern = paste0(data_to_process, "_subclustered.rds"),
    full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)
names(seurat_corticotrophs) <- gsub("(\\.rds|rds_outs/|_subclustered)", "", filenames)

pred_dir <- "label_transfer_model_output/predictions/"

for (dataset in unique(datasets$study_id)) {
    filename <- paste0(pred_dir, dataset, "_predictions.csv")
    predictions <- read.csv(filename)
    print(paste("Adding predictions to", dataset))
    # Add the predictions to the Seurat objects
    seurat_corticotrophs[[dataset]]@meta.data <- cbind(seurat_corticotrophs[[dataset]]@meta.data, predictions[3:ncol(predictions)])

    seurat_corticotrophs[[dataset]]@meta.data %>%
        mutate_at(vars(starts_with("Predicted")), as.factor) -> seurat_corticotrophs[[dataset]]@meta.data
}

confusion_matrices <- lapply(seurat_corticotrophs, function(obj) {
    predictions <- obj@meta.data %>% select(starts_with("Predicted"))
    predictors_names <- sub("Predicted_cluster_", "", names(predictions))

    # Get the confusion matrix for each dataset
    conf <- lapply(seq_along(predictions), function(i) {
        pred <- predictions[[i]]

        # Make also sure we have all of the classes of the predictor datasets
        if (!all(levels(Idents(seurat_corticotrophs[[predictors_names[i]]])) %in% levels(pred))) {
            pred <- factor(pred, levels = levels(Idents(seurat_corticotrophs[[predictors_names[i]]])))
        }

        # Make sure we have the -1 class (sometimes it's missing if all cells are predicted with high confidence)
        if (!-1 %in% levels(pred)) {
            pred <- factor(pred, levels = c(-1, levels(pred)))
        }

        tb <- table(Cluster1 = Idents(obj), Cluster2 = pred)
        # Convert to percentages
        # We normalise by column (i.e. by the predicted cluster)
        tb <- as.table(apply(tb, 2, function(x) x / sum(x)))
        tb
    })

    names(conf) <- predictors_names

    conf
})

# Convert the list of confusion matrices into a single data frame
confusion_matrices_df <- lapply(seq_along(confusion_matrices), function(i) {
    conf_df <- lapply(seq_along(confusion_matrices[[i]]), function(j) {
        predictor_name <- names(confusion_matrices[[i]])[j]

        conf <- confusion_matrices[[i]][[j]]
        print(conf)
        conf %>%
            as.data.frame() %>%
            mutate(
                Dataset1_name = names(confusion_matrices)[[i]],
                Dataset2_name = predictor_name
            ) -> conf_df
        conf_df
    })

    conf_df
})

# Note this is a list of lists, so we need to apply do.call twice
confusion_matrices_df <- do.call(rbind, lapply(confusion_matrices_df, do.call, what = rbind)) %>%
    filter(Cluster2 != -1) %>%
    filter(Dataset1_name != Dataset2_name) -> confusion_matrices_df

format_dataset_name <- function(name) {
    #' Go from Xxxxxxxxx-n to Xn
    #' @param name: The name to format
    #' @return: The formatted name
    #' @examples
    #' format_dataset_name("Lopez2021M-2") will return 'L2'
    #' format_dataset_name("Vennekens2021M-10") will return 'V10'
    paste0(substr(name, 1, 1), substr(name, nchar(name), nchar(name)))
}

# This is adapted from 06_compare_datasets.R
get_similarity_graph <- function(min_similarity, do_plot = TRUE,
                                 node_palette = NULL,
                                 label_color = "black",
                                 out_radius = 5, in_radius = 1) {
    #' Plots a graph, connecting subcluster with at least `min_similarity` cells predicted to be the same
    #' @param min_similarity: The minimum percentage of cells predicted to be the same
    #' @param do_plot: Whether to plot the graph - default is TRUE
    #' @param node_palette: The palette for the nodes - NULL for default (Set2)
    #' @param label_color: The colour for the labels - default is black

    # We find all of the combinations of dataset name and subcluster number
    datasets <- unique(c(
        paste(confusion_matrices_df$Dataset1_name,
            confusion_matrices_df$Cluster1,
            sep = "-"
        ),
        paste(confusion_matrices_df$Dataset2_name,
            confusion_matrices_df$Cluster2,
            sep = "-"
        )
    ))

    # Simplify the naming (for graphical purposes) and sort alphabetically
    datasets <- sort(format_dataset_name(datasets), decreasing = TRUE)
    edges <- confusion_matrices_df %>%
        filter(Freq >= min_similarity) %>%
        mutate(
            # Note the direction of the edges
            V1 = paste0(substr(Dataset2_name, 1, 1), Cluster2),
            V2 = paste0(substr(Dataset1_name, 1, 1), Cluster1)
        ) %>%
        select(V1, V2, Freq)

    # Now create a graph using the subcluster as vertices and connecting those
    # with >= % cells predicted to be from the same cluster
    # Note that this is a directed graph, contrary to the graph
    # in 06_compare_datasets.R (where we don't have direction
    # as the % of shared markers is symmetric)
    network <- graph_from_data_frame(
        d = edges,
        vertices = datasets,
        directed = TRUE
    )

    # Create a layout for our network

    data.frame(FullName = datasets) %>%
        mutate(
            Dataset = substr(FullName, 1, 1),
            SubCluster = substr(FullName, 2, 2)
        ) %>%
        group_split(Dataset) -> nodes_by_group

    # # Find communities
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
        weight <- E(network)$Freq - min_similarity
        plot(network,
            vertex.color = node_colors,
            edge.width = weight * 50,
            edge.curved = 0.3,
            edge.color = rgb(0, 0, 0, 0.3 * weight * (1 / (1 - min_similarity))),
            edge.arrow.size = 1.5,
            vertex.label.color = label_color,
            layout = layout,
            # layout = layout.fruchterman.reingold
            main = paste0(
                "Clusters with >= ", format(min_similarity * 100, digits = 2),
                "% similarity"
            )
        )
    }

    # Return the graph object
    return(network)
}

for (min_similarity in c(0.8, 0.9, 0.95)) {
    if (output_format == "png") {
        png(paste0(output_folder, "label_transfer_graph_", min_similarity, ".png"),
            width = 10, height = 10,
            units = "in", res = 300
        )
    } else if (output_format == "pdf") {
        pdf(paste0(output_folder, "label_transfer_graph_", min_similarity, ".pdf"),
            width = 10, height = 10
        )
    }

    community_palette <- c(
        "#F16745", "#FFC65D", "#7BC8A4",
        "#4CC3D9", "#93648D", "#808070"
    )

    graph <- get_similarity_graph(min_similarity,
        do_plot = TRUE,
        node_palette = community_palette
    )
    dev.off()

    memberships <- membership(walktrap.community(graph, steps = 100))

    tb <- table(memberships)

    # Assign a colour to each community.
    comm_colors <- community_palette[sort(unique(memberships))]
    # Communities with a single node are grey
    comm_colors[tb == 1] <- "lightgrey"
    names(comm_colors) <- names(tb)

    for (i in seq_along(seurat_corticotrophs)) {
        # Add a column to the metadata with the community
        seurat_corticotrophs[[i]]@meta.data[[paste0("nn_community_", min_similarity)]] <-
            memberships[match(
                paste0(
                    substr(seurat_corticotrophs[[i]]$author[1], 1, 1),
                    Idents(seurat_corticotrophs[[i]])
                ),
                names(memberships)
            )]
        seurat_corticotrophs[[i]]@meta.data[[paste0("nn_community_", min_similarity)]] <- factor(seurat_corticotrophs[[i]]@meta.data[[paste0("nn_community_", min_similarity)]], levels = sort(unique(memberships)))
    }

    # Now plot the UMAP reductions, colouring by community
    community_plots <- lapply(seurat_corticotrophs, function(s) {
        # Get the first letter of the dataset name.
        # This corresponds to the names of the graph nodes
        graph_nodes_name <- substr(s$author[1], 1, 1)

        palette <- comm_colors

        p <- DimPlot(s, group.by = paste0("nn_community_", min_similarity)) +
            scale_color_manual(
                name = "Community",
                values = palette
            ) +
            ggtitle(paste(s$author[1], s$year[1], "-", s$sex[1])) +
            theme(legend.position = "none")
    })

    if (output_format == "pdf") {
        pdf(paste0(output_folder, "umap_nn_communities_min_sim_", min_similarity, "_", data_to_process, ".pdf"),
            width = 10, height = 10
        )
    } else {
        png(paste0(output_folder, "umap_nn_communities_min_sim_", min_similarity, "_", data_to_process, ".png"),
            width = 10, height = 10, units = "in", res = 300
        )
    }

    grid.arrange(
        grobs = community_plots, ncol = 3,
        top = textGrob(
            paste("Cells with >= ", format(min_similarity * 100, digits = 2), "% similarity"),
            gp = gpar(fontsize = 20)
        )
    )
    dev.off()
} # end for min_similarity
