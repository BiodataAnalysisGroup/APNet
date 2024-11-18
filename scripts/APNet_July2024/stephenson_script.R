library(Seurat)
setwd("/Users/georgegavriilidis/pertpy/output_folder")

Raw_data <- Seurat::Read10X(data.dir = 'matrix_files')
metadata <- read.csv('metadata.csv')
sobj <- CreateSeuratObject(counts = Raw_data, meta.data = metadata)
embedR <- read.csv('embedding_umap.csv')
embedR <- as.matrix(embedR)
rownames(embedR) <- embedR[, 1]  # Use the first column as row names
embedR <- embedR[, -1]
mode(embedR) <- "numeric"
sobj[["umap"]] <- CreateDimReducObject(embedR, key = "umap_", assay = "RNA")
DimPlot(sobj, group.by = "Status")
saveRDS(sobj, "stephenson_seurat_covid19.rds")
