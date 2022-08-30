library(Seurat)
library(pbapply)
library(R.utils)

filenames <- dir("rds_outs", pattern = "subclustered.rds", 
                 full.names = TRUE)

seurat_obj_cort <- pblapply(filenames, readRDS)

# Find common genes
common_genes <- NULL

for (obj in seurat_obj_cort)
  {
  if (!length(common_genes))
    {
    common_genes <- rownames(obj)
    }
  else
    {
    common_genes <- intersect(common_genes, rownames(obj))
    }
  }

if (!dir.exists("expr_matrices"))
  dir.create("expr_matrices")

pbsapply(seurat_obj_cort, function(obj){
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], ".csv")

  # Only save genes that are common to all datasets
  write.csv(GetAssayData(obj)[common_genes,], file = outfile,
         row.names = TRUE)
  gzip(outfile, overwrite = TRUE)
  
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_clusters.csv")
  write.csv(data.frame(Barcode = Cells(obj), Cluster = Idents(obj)), 
            outfile, row.names = FALSE, quote = FALSE)
  })

