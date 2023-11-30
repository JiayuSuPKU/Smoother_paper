library(Seurat)
library(SeuratDisk)
library(tidyverse)

# load single-cell reference data
CRC_SingleCell <- readRDS("~/Projects/Smoother_paper/data/crc_stereo/CRC_SingleCell.rds")

# save raw counts and convert into h5ad
scref_raw = CreateSeuratObject(
  counts = GetAssayData(CRC_SingleCell, 'counts'),
  meta.data = data.frame('orig.ident' = CRC_SingleCell@meta.data$orig.ident,
                         'cell_type_orig' = as.character(Idents(CRC_SingleCell)),
                         row.names = rownames(CRC_SingleCell@meta.data))
)
SaveH5Seurat(scref_raw, filename = "~/Projects/Smoother_paper/data/crc_stereo/scref/scref_raw.h5Seurat", overwrite = T)
Convert("~/Projects/Smoother_paper/data/crc_stereo/scref/scref_raw.h5Seurat", dest = "h5ad", overwrite = T)

# load stereo-seq spatial data
CRC_Stereo <- readRDS("~/Projects/Smoother_paper/data/crc_stereo/CRC_Stereo.rds")

# save stereo-seq data of sample P19_T (tumor) and convert into h5ad
p19_t <- CreateSeuratObject(
  counts = GetAssayData(subset(CRC_Stereo, subset = id == 'P19_T'))
)
p19_t <- AddMetaData(p19_t, 
                     subset(CRC_Stereo, subset = id == 'P19_T')@meta.data$bayes_clusters,
                     'bayes_clusters')
coords <- subset(CRC_Stereo, subset = id == 'P19_T')@images$image_P19_T@coordinates[,c('x', 'y')]
p19_t <- AddMetaData(p19_t, coords)
SaveH5Seurat(p19_t, filename = "~/Projects/Smoother_paper/data/crc_stereo/P19_T.h5Seurat", overwrite = T)
Convert("~/Projects/Smoother_paper/data/crc_stereo/P19_T.h5Seurat", dest = "h5ad", overwrite = T)

# save stereo-seq data of sample P19_NT (normal tissue) and convert into h5ad
p19_nt <- CreateSeuratObject(
  counts = GetAssayData(subset(CRC_Stereo, subset = id == 'P19_NT'))
)
p19_nt <- AddMetaData(p19_nt, 
                     subset(CRC_Stereo, subset = id == 'P19_NT')@meta.data$bayes_clusters,
                     'bayes_clusters')
coords <- subset(CRC_Stereo, subset = id == 'P19_NT')@images$image@coordinates[,c('x', 'y')]
p19_nt <- AddMetaData(p19_nt, coords)
SaveH5Seurat(p19_nt, filename = "~/Projects/Smoother_paper/data/crc_stereo/P19_NT.h5Seurat", overwrite = T)
Convert("~/Projects/Smoother_paper/data/crc_stereo/P19_NT.h5Seurat", dest = "h5ad", overwrite = T)
