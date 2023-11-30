library(Seurat)
library(SeuratDisk)
library(tidyverse)

data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/prostate_ref_mapping/"

# store healthy ST slides
slide.seq.raw.counts <- readRDS(paste0(data_dir, "slide.seq.raw.counts.rds"))
slide.seq.ano <- readRDS(paste0(data_dir, "slide.seq.ano.rds"))
slide.seq.BeadLocations <- readRDS(paste0(data_dir, "slide.seq.BeadLocations.rds"))

for (sample in c('HP1', 'HP2', 'HP3', 'HP4')){
  st.mtx <- slide.seq.raw.counts[[sample]]
  st.anno <- cbind(
    slide.seq.ano[colnames(st.mtx),], 
    slide.seq.BeadLocations[[sample]]
  ) %>% mutate(batch = paste0(sample, '_ST'))
  
  st <- CreateSeuratObject(
    counts = st.mtx,
    meta.data = st.anno,
    min.features = 10
  )

  fname = sprintf("%s/ST_%s.h5Seurat", data_dir, sample)
  SaveH5Seurat(st, filename = fname, overwrite = T)
  Convert(fname, dest = "h5ad", overwrite = T)
}


