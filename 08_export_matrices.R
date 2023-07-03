library(Seurat)
library(Matrix)
library(R.utils)

filenames <- dir("rds_outs",
  pattern = "subclustered.rds",
  full.names = TRUE
)

seurat_obj_cort <- pblapply(filenames, readRDS)

if (!dir.exists("expr_matrices")) {
  dir.create("expr_matrices")
}

sapply(seurat_obj_cort, function(obj) {
  print(paste(obj$author[1], obj$year[1], "-", obj$sex[1]))
  print("---------------------------------")

  print("Exporting expression matrix")
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_expression.csv")
  write.csv(GetAssayData(obj), file = outfile, row.names = TRUE)
  # gzip(outfile, overwrite = TRUE)

  print("Exporting raw counts")
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_counts")
  counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
  write.csv(counts_matrix, file = paste0(outfile, ".csv"), row.names = TRUE)
  writeMM(counts_matrix, file = paste0(outfile, ".mtx"))

  print("Exporting cluster assignments")
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_clusters.csv")
  write.csv(
    data.frame(
      Barcode = Cells(obj),
      Cluster = Idents(obj)
    ),
    outfile,
    row.names = FALSE, quote = FALSE
  )

  print("Exporting dimension reductions (PCA and UMAP)")
  outfile_pca <- paste0("expr_matrices/", obj$orig.ident[1], "_pca.csv")
  outfile_umap <- paste0("expr_matrices/", obj$orig.ident[1], "_umap.csv")
  pca <- Embeddings(obj, reduction = "pca")
  umap <- Embeddings(obj, reduction = "umap")
  write.csv(cbind(Barcode = Cells(obj), pca),
    outfile_pca,
    row.names = FALSE, quote = FALSE
  )
  write.csv(data.frame(Barcode = Cells(obj), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2]),
    outfile_umap,
    row.names = FALSE, quote = FALSE
  )

  print("Exporting metadata")
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_metadata.csv")
  write.csv(obj@meta.data, outfile, row.names = FALSE, quote = FALSE)
})

# # Get all SCT filenames
# filenames <- dir("rds_outs", pattern = "M_SCT.rds", full.names = TRUE)
# # Remove Ho 2020
# filenames <- filenames[!grepl("Ho2020", filenames)]

# print("Exporting expression matrices for other cells")

# pbsapply(filenames, function(filename) {
#   # Now save a matrix containing non-corticotrophs pituitary cells.
#   # This is used when determining softmax threshold
#   obj <- readRDS(filename)
#   filename_cort <- gsub("_SCT.rds", "_corticotrophs.rds", filename)
#   obj_cort <- readRDS(filename_cort)

#   all_cells <- Cells(obj)
#   other_cells <- setdiff(all_cells, Cells(obj_cort))

#   obj <- subset(obj, cells = other_cells)
#   outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_other.csv")

#   # Check we have all the common genes
#   if (!all(common_genes %in% rownames(obj))) {
#     stop(paste("Missing genes in file", filename))
#   }

#   write.csv(GetAssayData(obj)[common_genes, ], outfile, row.names = TRUE)
#   gzip(outfile, overwrite = TRUE)
# })
