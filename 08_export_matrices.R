library(Seurat)
library(pbapply)
library(R.utils)

filenames <- dir("rds_outs",
  pattern = "subclustered.rds",
  full.names = TRUE
)

seurat_obj_cort <- pblapply(filenames, readRDS)

# Find common genes
common_genes <- NULL

for (obj in seurat_obj_cort)
{
  if (!length(common_genes)) {
    common_genes <- rownames(obj)
  } else {
    common_genes <- intersect(common_genes, rownames(obj))
  }
}

if (!dir.exists("expr_matrices")) {
  dir.create("expr_matrices")
}

print("Exporting corticotrophs expression matrices")

pbsapply(seurat_obj_cort, function(obj) {
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], ".csv")

  # Only save genes that are common to all datasets
  write.csv(GetAssayData(obj)[common_genes, ],
    file = outfile, row.names = TRUE
  )
  gzip(outfile, overwrite = TRUE)

  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_clusters.csv")
  write.csv(data.frame(Barcode = Cells(obj), Cluster = Idents(obj)),
    outfile,
    row.names = FALSE, quote = FALSE
  )
})

# Get all SCT filenames
filenames <- dir("rds_outs", pattern = "M_SCT.rds", full.names = TRUE)
# Remove Ho 2020
filenames <- filenames[!grepl("Ho2020", filenames)]

print("Exporting expression matrices for other cells")

pbsapply(filenames, function(filename) {
  # Now save a matrix containing non-corticotrophs pituitary cells.
  # This is used when determining softmax threshold
  obj <- readRDS(filename)
  filename_cort <- gsub("_SCT.rds", "_corticotrophs.rds", filename)
  obj_cort <- readRDS(filename_cort)

  all_cells <- Cells(obj)
  other_cells <- setdiff(all_cells, Cells(obj_cort))

  obj <- subset(obj, cells = other_cells)
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_other.csv")

  # Check we have all the common genes
  if (!all(common_genes %in% rownames(obj))) {
    stop(paste("Missing genes in file", filename))
  }

  write.csv(GetAssayData(obj)[common_genes, ], outfile, row.names = TRUE)
  gzip(outfile, overwrite = TRUE)
})