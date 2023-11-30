
library(ArchR)
library(Seurat)
library(grid)
library(Signac)
library(parallel)
library(grid)
library(ggplot2)
library(patchwork)
library(dplyr)
source('process_code/Data_visualization/SpatialDimPlot_new.R')
source('process_code/Data_visualization/getGeneScore_ArchR.R')
source('process_code/Data_visualization/SpatialPlot_new.R')
threads = 20
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- './GSM5028436_ME11_H3K27ac_50um.fragments.tsv.gz'
sampleNames <- 'ME11_H3K27ac_50um'

## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj


## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data.dir, "GSM5028436_spatial"), 
                       filter.matrix = filter.matrix)
meta.data.spatial <- meta.data[paste('ME11_H3K27ac_50um#', row.names(image@coordinates), '-1', sep=""), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue

proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)


proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)

proj_in_tissue <- addImputeWeights(proj_in_tissue)

markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList_neg <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC <= -1")

# Export marker genes
for (i in seq_len(length(markerList_pos))) {
  write.table(markerList_pos[[i]], file=paste0('./markers_list/', sampleNames, '_C', i, '_markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)


proj_in_tissue <- addGroupCoverages(proj_in_tissue)
proj_in_tissue <- addReproduciblePeakSet(proj_in_tissue)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)
meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names


gene_score <- getGeneScore_ArchR(colorBy = "GeneScoreMatrix", ArchRProj = proj_in_tissue, name = markerGenes, imputeWeights = getImputeWeights(proj_in_tissue))

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)

image <- Read10X_Image(image.dir = file.path(data.dir, "GSM5028436_spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

objectooio
''

save.image("H3K27ac.50um.RData")

load("H3K27ac.50um.RData")
refrna = read.csv("processed/ref.var_genes.bycluster.csv")$X
gmat = getMatrixFromProject(proj_in_tissue)
refrna = refrna[refrna %in% gmat@elementMetadata$name]
# gene_score <- getGeneScore_ArchR(colorBy = "GeneScoreMatrix", ArchRProj = proj_in_tissue, name = refrna, imputeWeights = getImputeWeights(proj_in_tissue))
all_gene_score <- getGeneScore_ArchR(colorBy = "GeneScoreMatrix", ArchRProj = proj_in_tissue, name =  gmat@elementMetadata$name, imputeWeights = getImputeWeights(proj_in_tissue))
write.csv(all_gene_score, file = "processed/H3K27ac.50um.all_gene_score.csv")
var_gene_score <- all_gene_score[order(rowVars(as.matrix(all_gene_score)),decreasing = T)[1:500],]
write.csv(var_gene_score, file = "processed/H3K27ac.50um.var_gene_score.csv")
# write.csv(gene_score, file = "processed/H3K27ac.50um.gene_score.csv")
write.csv(image@coordinates[,4:5], file = "processed/H3K27ac.50um.coord.csv")

