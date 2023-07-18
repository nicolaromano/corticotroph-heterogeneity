#' 11_prepare_velocity_analysis.R'
#' This script prepares data for velocity analysis using scVelo
#' It goes through all of the samples and creates a barcodes.tsv file
#' containing only the barcodes for corticotrophs. This will speed up the calculations substantially.

library(Seurat)
library(pbapply)
library(velocyto.R)
library(dplyr)

datasets <- read.csv("datasets.csv")

filenames <- dir("rds_outs",
    pattern = "subclustered.rds",
    full.names = TRUE
)
seurat_objects <- pblapply(filenames, readRDS)

names(seurat_objects) <- sapply(seurat_objects, function(obj) {
    return(obj$orig.ident[1])
})

sapply(1:nrow(datasets), function(row) {
    sample <- datasets[row, ]
    barcodes_file <- paste0(sample["cr_output"], "/barcodes.tsv.gz")
    all_barcodes <- read.csv(barcodes_file, header = FALSE)
    study <- unlist(sample["study_id"])
    sample <- unlist(sample["source"])
    print(paste0("Processing ", study, "::", sample))

    cort_barcodes <- Cells(seurat_objects[[study]])

    # print(head(cort_barcodes))
    # print(head(all_barcodes))
    filtered_barcodes <- intersect(cort_barcodes, unlist(all_barcodes))
    print(paste("Out of", length(all_barcodes$V1), "barcodes,", length(filtered_barcodes), "are from corticotrophs"))
    print(head(filtered_barcodes))
    outfile <- sub(barcodes_file, pattern = "barcodes.tsv.gz", replacement = "cort_barcodes.tsv")
    write.table(filtered_barcodes, file = outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
})

stop("RUN VELOCYTO NOW!")

# Velocyto should be run now, to generate the loom files for example
# velocyto run --bcfile <samplepath>/filtered_feature_bc_matrix/cort_barcodes.tsv \
#     --sampleid <samplename> \
#     --outputfolder RNA_velocity/<samplename>/ \
#     --samtools-threads 16 --samtools-memory 4096 \
#     <samplepath>/possorted_genome_bam.bam <genomepath>/gencode.vM23.primary_assembly.annotation.gtf.filtered

# At this point, we could process the loom files if we wanted to do velocity analysis using velocyto.R
# Currently this simply crashes R without any error message, so we'll use scVelo instead.
# To use scVelo, use the 12_scvelo.ipybn notebook.

# # From https://stackoverflow.com/questions/46079152/how-to-merge-big-sparse-matrices
# merge_sparse <- function(M.list) {
#     #' Merges a list of sparse matrices
#     #' @param M list of sparse matrices
#     #' @return a sparse matrix

#     A <- M.list[[1]]

#     for (i in 2:length(M.list)) { # i indexes of matrices
#         # finding what's missing
#         misA <- colnames(M.list[[i]])[!colnames(M.list[[i]]) %in% colnames(A)]
#         misB <- colnames(A)[!colnames(A) %in% colnames(M.list[[i]])]

#         misAl <- as.vector(numeric(length(misA)), "list")
#         names(misAl) <- misA
#         misBl <- as.vector(numeric(length(misB)), "list")
#         names(misBl) <- misB

#         ## adding missing columns to initial matrices
#         An <- Reduce(cbind, c(A, misAl))
#         if (length(misA) > 0) {
#             lenA <- ncol(An) - length(misA) + 1
#             colnames(An)[lenA:ncol(An)] <- names(misAl)
#         }

#         Bn <- Reduce(cbind, c(M.list[[i]], misBl))
#         if (length(misB) > 0) {
#             lenB <- ncol(Bn) - length(misB) + 1
#             colnames(Bn)[lenB:ncol(Bn)] <- names(misBl)
#         }

#         Bn <- Bn[, colnames(An)]

#         # final bind
#         A <- rbind(An, Bn, use.names = T)
#         print(paste(i, "/", length(M.list)))
#     }
#     A
# }

# process_datasets <- function(datasets, loom_dir, skip_existing = FALSE) {
#     #' Processes the loom files for each dataset
#     #'
#     #' @param datasets A data.frame containing the study_id and source columns
#     #' @param loom_dir The directory containing the loom files
#     #' @param skip_existing Whether to skip datasets for which the output files already exist
#     #'
#     #' @return Nothing, but writes the output files to disk as RDS files

#     sapply(unique(datasets$study_id), function(dataset) {
#         print(paste0("Processing ", dataset))
#         print("--------------------")

#         # Check all loom files are present
#         sources <- datasets %>%
#             filter(study_id == dataset) %>%
#             select(study_id, source)

#         filenames <- paste0(loom_dir, "/", sources$study_id, "_", sources$source, ".loom")
#         if (any(!file.exists(filenames))) {
#             print("The following files are missing:")
#             # Print missing files
#             print(filenames[!file.exists(filenames)])
#             print("Skipping this dataset")
#             return(NULL)
#         }

#         # Check if output files already exist
#         spliced_file <- paste0(loom_dir, "/", dataset, "_spliced.rds")
#         unspliced_file <- paste0(loom_dir, "/", dataset, "_unspliced.rds")
#         amb_file <- paste0(loom_dir, "/", dataset, "_ambiguous.rds")
#         if (skip_existing && file.exists(spliced_file) && file.exists(unspliced_file) && file.exists(amb_file)) {
#             print(paste("Output files already exist. Skipping dataset", dataset))
#             return(NULL)
#         }

#         mat <- pblapply(filenames, function(filename) {
#             return(read.loom.matrices(filename))
#         })

#         spliced <- lapply(mat, function(x) {
#             return(x$spliced)
#         })

#         unspliced <- lapply(mat, function(x) {
#             return(x$unspliced)
#         })

#         ambiguous <- lapply(mat, function(x) {
#             return(x$ambiguous)
#         })

#         if (length(spliced) > 1) {
#             print("Merging sparse matrices...")
#             all_spliced <- merge_sparse(spliced)
#             all_unspliced <- merge_sparse(unspliced)
#             all_amb <- merge_sparse(ambiguous)
#         } else {
#             all_spliced <- spliced[[1]]
#             all_unspliced <- unspliced[[1]]
#             all_amb <- ambiguous[[1]]
#         }

#         saveRDS(all_spliced, spliced_file)
#         saveRDS(all_unspliced, unspliced_file)
#         saveRDS(all_amb, amb_file)
#     })

#     return(invisible(NULL))
# }

# process_datasets(datasets, "RNA_velocity", skip_existing = TRUE)