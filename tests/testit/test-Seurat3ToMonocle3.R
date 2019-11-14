library(testit)
library(Seurat)

assert('Seurat3ToMoncle3 returns a cell_', {
  pbmc_raw <- read.table(file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),as.is = TRUE)
  pbmc_small <- CreateSeuratObject(counts = pbmc_raw)
  pbmc_monocle3 <- Seurat3ToMonocle3(pbmc_small)
  as.character(class(pbmc_monocle3)) %==% "cell_data_set"

})
