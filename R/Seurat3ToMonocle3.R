##' @title convert different formats from different softwares or packages
##' @description This is an package to convert formats among different softwares or packages
##'
##' @details nothing
##' @param seuObj the input seurat3 object
##' @return an cell_data_set object
##' @importFrom monocle3 new_cell_data_set
##' @export Seurat3ToMonocle3
##' @examples
##' pbmc_raw <- read.table(file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),as.is = TRUE)
##' pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
##' pbmc_monocle3 <- Seurat3ToMonocle3(pbmc_small)
##'
##' @author Feng Zhang

Seurat3ToMonocle3 = function(seuObj){
  cds <- new_cell_data_set(as(seuObj@assays$RNA@data, "sparseMatrix"),
                           cell_metadata = as.data.frame(seuObj@meta.data),
                           gene_metadata = data.frame(gene_short_name = rownames(seuObj), row.names = rownames(seuObj))
  )
  #assign cluster information from seurat object
  recreate_partition <- factor(rep(NA, length(cds@colData@rownames)))
  names(recreate_partition) <- cds@colData@rownames
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate_partition

  if(!is.null(seuObj@meta.data$seurat_clusters)) {
    seurat_clusters=factor(seuObj@meta.data$seurat_clusters)
    names(seurat_clusters) = row.names(seuObj@meta.data)
    cds@clusters@listData[["UMAP"]][["clusters"]] = seurat_clusters
  }


  # assign reduced dimmension
  if(!is.null(seuObj@reductions[["umap"]])){
    cds@reducedDims@listData[["UMAP"]] <- seuObj@reductions[["umap"]]@cell.embeddings
  }

  if(!is.null(seuObj@reductions[["pca"]])){
    cds@reducedDims@listData[["PCA"]] <- seuObj@reductions[["pca"]]@cell.embeddings
  }

  cds
}
