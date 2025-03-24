library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pbapply)
library(KEGGREST)
library(biomaRt)
library(boot)

data_to_process <- "F"
datasets <- read.csv("datasets.csv")

output_format <- "pdf" # "pdf", "png" or "none"

filenames <- dir("rds_outs",
    pattern = paste0(data_to_process, "_subclustered.rds"),
    full.names = TRUE
)

seurat_corticotrophs <- pblapply(filenames, readRDS)

names(seurat_corticotrophs) <- sapply(seurat_corticotrophs, function(s) {
    s$orig.ident[1]
})

get_KEGG_pathway_genes <- function(pathway, cache = TRUE) {
    #' Get genes from KEGG pathway
    #' @param pathway KEGG pathway ID
    #' @param cache If TRUE, save the result to a file
    #' @return Data frame with columns ID, Gene and Description

    if (cache) {
        # Check if the pathway has already been downloaded
        if (file.exists(paste0("kegg_", pathway, ".rds"))) {
            p <- readRDS(paste0("kegg_", pathway, ".rds"))
        } else {
            p <- keggGet(pathway)
            saveRDS(p, paste0("kegg_", pathway, ".rds"))
        }
    }

    # p[[1]]$GENE is a character vector with entries
    # GeneID, GeneName;GeneDescription, GeneID, GeneName;GeneDescription, ...
    # We first split the vector into a matrix with two columns, so that we can
    # isolate the ID from the rest of the information
    # We then split the second column into two columns, so that we can separate
    # gene name and description
    res <- matrix(p[[1]]$GENE, ncol = 2, byrow = TRUE) %>%
        as.data.frame() %>%
        setNames(c("ID", "Gene")) %>%
        # fill = "left" is needed to handle cases where the gene name is missing
        separate(Gene, into = c("Gene", "Description"), sep = ";", fill = "left") %>%
        filter(!is.na(Gene))

    return(res)
}

get_GO_genes <- function(go_term) {
    # Get genes from GO term
    ensembl <- useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "mmusculus_gene_ensembl"
    )

    go_term <- getBM(
        attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
        filters = "go", values = go_term, mart = ensembl
    )
    colnames(go_term) <- c("ID", "Gene")
    return(go_term)
}

# To find specific pathways you can use e.g.
# grep("immediate", keggList("pathway", "mmu"), value = TRUE, ignore.case = TRUE)

plot_pathway_scores <- function(pathways, extra_gene_list = NULL, title = "") {
    #' Plot module scores for a given pathway
    #' @param pathways list of KEGG pathway ID or GO term.    
    #' @param extra_gene_list Extra list of genes to include
    #' @param title Title of the plot
    #' @return ggplot object

    genes <- c()

    for (p in pathways) {
        if (grepl("^GO", p)) {
            genes <- c(genes, get_GO_genes(p)$Gene)
        } else {
            genes <- c(genes, get_KEGG_pathway_genes(p)$Gene)
        }
    }

    if (length(extra_gene_list) > 0) {
        genes <- c(genes, extra_gene_list)
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

    # Only keep genes that are present in the Seurat object
    genes <- genes[genes %in% rownames(seurat_corticotrophs[[1]])]
    if (length(genes) == 0) {
        stop("No genes found in Seurat object")
    }
    # Add module score to each object
    objects_with_score <- lapply(seurat_corticotrophs, function(obj) {
        AddModuleScore(obj, features = list(genes), name = "pathway", nbin = 10)
    })

    df <- lapply(seq_along(objects_with_score), function(i) {
        obj <- objects_with_score[[i]]
        return(data.frame(
            Community = as.numeric(obj$marker_community_10),
            Pathway = obj[["pathway1"]],
            Dataset = paste(obj$author, obj$year)
        ))
    })

    df <- do.call(rbind, df)
    # Note that when we extract the module score we get a named vector, which overwrite the column name
    # We need to rename the columns for the following code to work
    colnames(df) <- c("Community", "Pathway", "Dataset")

    # Remove NA communities
    df <- df %>% filter(!is.na(Community))
    df$Community <- factor(df$Community)

    g <- ggplot(df, aes(x = Community, y = Pathway, group = Community)) +
        geom_violin(aes(fill = Community)) +
        scale_fill_manual(values = community_palette, name = "Community") +
        stat_summary(fun = median, geom = "point", color = "black", size = 3) +
        xlab("Community") +
        ylab("Pathway score") +
        ggtitle(title) +
        theme_minimal() +
        theme(
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        )

    return(list(plot = g, data = df, genes = genes))
}

# List from Tullai et al. JBC 2007
IEG <- c(
    "Atf3", "Bhlhe40", "Ccl2", "Ccn1", "Ccn2", "Ccnl1", "Cebpd",
    "Csrnp1", "Cxcl3", "Dusp1", "Dusp5", "Dusp6", "Egr1", "Egr3",
    "F3", "Flg", "Fos", "Fosb", "Gadd45b", "Gbp1", "Gem", "Hbegf",
    "Ier3", "Il6", "Jun", "Junb", "Klf10", "Klf6", "Klhl21", "Ldlr",
    "Mcl1", "Nfkbia", "Nfkbiz", "Nr2c2", "Nr4a1", "Nr4a2", "Plau",
    "Pmaip1", "Rcan1", "Sgk1", "Slc2a3", "Srf", "Tnfaip3", "Trib1", "Tsc22d1", "Zfp36"
)

p1 <- plot_pathway_scores(c("GO:0009725", "mmu04081"), title = "Response to hormone")
p2 <- plot_pathway_scores(NULL, extra_gene_list = IEG, title = "Immediate early genes")
p3 <- plot_pathway_scores("mmu00190", title = "Oxidative phosphorylation")
p4 <- plot_pathway_scores("GO:0006099", title = "TCA cycle")


if (output_format == "pdf") {
    pdf(paste0(
        "plots/pathway_scores_markers_communities", data_to_process,
        ".pdf"
    ), width = 12, height = 8)
} else if (output_format == "png") {
    png(paste0(
        "plots/pathway_scores_markers_communities", data_to_process,
        ".png"
    ), width = 12, height = 8, units = "in", res = 300)
}

p1$plot + p2$plot + p3$plot + p4$plot

if (output_format != "none") {
    dev.off()
}

######## STATS ########


# Analysing this data can be quite tricky because:
# - The data is not normally distributed and often has a multimodal distribution, but not in all groups
# - The data is not homoscedastic
# - The data is not independent, so we need to account for the nested structure of the data (cells within datasets) to avoid pseudoreplication
# - The data is not balanced

# We use a bootstrap approach to estimate the p-value for the difference between each pair of communities
# This is more robust than using a parametric test (e.g. mixed model) because it does not rely on the assumptions of the model

boot_median_diff <- function(data, indices) {
    #' Calculate the difference in median between two groups, stratified by dataset
    #' @param data Data frame with columns Community, Pathway, Dataset
    #' @param indices Indices of the data to use (resampled per dataset)
    #' @return Difference in median

    d <- data[indices, ]

    # Resample within each dataset
    d_resampled <- d %>%
        group_by(Dataset) %>%
        sample_n(size = n(), replace = TRUE) %>%
        ungroup()

    medians <- d_resampled %>%
        group_by(Community) %>%
        summarise(median = median(Pathway), .groups = "drop") %>%
        arrange(Community) %>%
        pull(median)

    return(medians[1] - medians[2])
}

get_boot_comparison <- function(df, nboot = 1000) {
    #' Get bootstrapped comparison between all pairs of communities
    #' Uses bias-corrected accelerated (BCa) confidence intervals
    #' @param df Data frame with columns Community and Pathway
    #' @param nboot Number of bootstrap samples
    #' @return Data frame with columns Community1, Community2, CI_low, CI_high, p_value

    communities <- as.numeric(as.character(sort(unique(df$Community))))
    comparisons <- expand.grid(Community1 = communities, Community2 = communities) %>%
        filter(Community1 < Community2)

    res <- lapply(1:nrow(comparisons), function(i) {
        df_filtered <- df %>%
            filter(Community %in% c(
                comparisons$Community1[i],
                comparisons$Community2[i]
            ))

        boot_res <- boot(
            df_filtered,
            statistic = boot_median_diff,
            R = nboot,
            parallel = "multicore", ncpus = parallel::detectCores()
        )

        CI <- boot.ci(boot_res, type = "bca")$bca[4:5]
        CI_low <- CI[1]
        CI_high <- CI[2]
        p_value <- 2 * min(
            mean(boot_res$t > 0, na.rm = TRUE),
            mean(boot_res$t < 0, na.rm = TRUE))

        return(data.frame(
            Community1 = comparisons$Community1[i],
            Community2 = comparisons$Community2[i],
            CI_low = CI_low,
            CI_high = CI_high,
            p_value = p_value)
        )
    })

    res <- do.call(rbind, res)
    return(res)
}

r1 <- get_boot_comparison(p1$data, 1000)
r2 <- get_boot_comparison(p2$data, 1000)
r3 <- get_boot_comparison(p3$data, 1000)
r4 <- get_boot_comparison(p4$data, 1000)

r <- rbind(r1, r2, r3, r4)

write.csv(r, paste0("tables/pathway_scores_markers_communities", data_to_process, ".csv"), row.names = FALSE)

