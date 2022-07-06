library(Seurat)
library(pbapply)
library(R.utils)

filenames <- dir("rds_outs", pattern = "subclustered.rds", 
                 full.names = TRUE)

seurat_obj_cort <- pblapply(filenames, readRDS)

if (!dir.exists("expr_matrices"))
  dir.create("expr_matrices")

pbsapply(seurat_obj_cort, function(obj){
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], ".csv")

  write.csv(GetAssayData(obj), file = outfile,
         row.names = TRUE)
  gzip(outfile, overwrite = TRUE)
  
  outfile <- paste0("expr_matrices/", obj$orig.ident[1], "_clusters.csv")
  write.csv(data.frame(Barcode = Cells(obj), Cluster = Idents(obj)), 
            outfile, row.names = FALSE, quote = FALSE)
  })
