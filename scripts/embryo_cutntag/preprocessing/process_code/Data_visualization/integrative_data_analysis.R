library(ArchR)
library(Seurat)
library(grid)

source('process_code/Data_visualization/SpatialDimPlot_new.R')

### Integrate with scRNA-seq reference data 
## https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
## Prepare MOCA data
MOCA_dir <- "./ref_data/MOCA/"

meta.data.RNA <- read.csv(file = paste0(MOCA_dir, 'cell_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- read.csv(file = paste0(MOCA_dir, 'gene_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- gene.ANN.RNA[, 'gene_short_name', drop = FALSE]

cds <- readRDS(paste0(MOCA_dir, 'gene_count_cleaned_sampled_100k.RDS'))

MOCA <- CreateSeuratObject(counts = cds, project = 'MOCA')
meta.data.RNA <- meta.data.RNA[colnames(MOCA), ]
meta.data.RNA <- meta.data.RNA[, c('Main_cell_type', 'development_stage')]

MOCA <- AddMetaData(object = MOCA, metadata = meta.data.RNA)
MOCA_E11 <- subset(MOCA, development_stage == 11.5)
MOCA_E11.raw.data <- as.matrix(GetAssayData(MOCA_E11, slot = 'data'))
MOCA_E11.raw.data <- as.data.frame(MOCA_E11.raw.data)
MOCA_E11.raw.data <- merge(gene.ANN.RNA, MOCA_E11.raw.data, by=0, all=TRUE)
which(is.na(MOCA_E11.raw.data$gene_short_name))

tt <- table(MOCA_E11.raw.data$gene_short_name)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x){
  rows <- which(MOCA_E11.raw.data$gene_short_name == x)
  return(rows[2:length(rows)] )
}
row_del <- unlist(lapply(name_rep, row_del_fun))
MOCA_E11.raw.data <- MOCA_E11.raw.data[-row_del, ]

row.names(MOCA_E11.raw.data) <- MOCA_E11.raw.data$gene_short_name
MOCA_E11.raw.data <- MOCA_E11.raw.data[, -c(1:2), drop=FALSE]
MOCA_E11 <- CreateSeuratObject(counts = MOCA_E11.raw.data, project = 'MOCA_E11', meta.data = MOCA_E11@meta.data)

MOCA_E11 <- FindVariableFeatures(MOCA_E11, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(MOCA_E11)
MOCA_E11 <- ScaleData(MOCA_E11, features = all.genes)
MOCA_E11 <- RunPCA(MOCA_E11, features = VariableFeatures(object = MOCA_E11))

MOCA_E11 <- FindNeighbors(MOCA_E11, dims = 1:10)
MOCA_E11 <- FindClusters(MOCA_E11, resolution = 0.1)
head(Idents(MOCA_E11), 5)
MOCA_E11 <- RunUMAP(MOCA_E11, dims = 1:10)

pdf("MOCA_E11.umap.pdf", 5, 5)
DimPlot(MOCA_E11, reduction = "umap", group.by = "Main_cell_type") + NoLegend()
dev.off()


write.csv((GetAssayData(object, 'data')[VariableFeatures(object = object),]), file='processed/20um.h3k27ac.var_genes.csv')

write.csv(AverageExpression(MOCA_E11, assays='RNA', group.by = 'Main_cell_type', slot='count')$RNA[VariableFeatures(object = MOCA_E11),], file='processed/ref.var_genes.bycluster.csv')



## Integration with ArchR oject
proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = MOCA_E11,
  addToArrow = TRUE,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)


## Plot results
meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('predictedCell', 'predictedGroup', 'predictedScore')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

Idents(spatial.obj) <- 'predictedGroup'

ids.highlight <- names(table(spatial.obj$predictedGroup))[3]
ids.highlight

p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = ids.highlight), 
                        facet.highlight = TRUE, pt.size.factor = 2.5, alpha = c(1,0), stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p


bin_counts = getMatrixFromProject(
  proj_in_tissue,
  binarize = TRUE,
  useMatrix = "TileMatrix")
rownames(bin_counts) = paste(bin_counts@elementMetadata$seqnames,as.character(bin_counts@elementMetadata$start+1),as.character(bin_counts@elementMetadata$start+5000), sep='-')
### Integrate with scCUT&Tag reference data (H3K4me3)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163532
K4.spatial <- CreateSeuratObject(counts = assays(bin_counts)$TileMatrix, assay = assay, meta.data = meta.data)

# K4.spatial      <- object #readRDS(file='/data/proj/GCB_MB/spatial_cut-tag/results/CT/H3K4me3/clustering/Seurat_object.Rds')
K4.spatial      <- RenameAssays(object = K4.spatial,bin_5000='bins_5000')


K4.spatial <- FindTopFeatures(K4.spatial, min.cutoff = "q0")
K4.spatial <- RunSVD(K4.spatial)
K4.spatial <- RunUMAP(K4.spatial, reduction = "lsi", dims = 2:30, reduction.name = "umap.spacial", reduction.key = "spatialUMAP_")

p1 <- DimPlot(K4.spatial)
p2 <- DimPlot(K4.spatial,group.by = 'YD.clusters')
p1+p2

K4.single_cell <- readRDS(file='ref_data/CT/h3k27ac.rds')

assay = 'bins_5000'

DefaultAssay(K4.spatial)      <- assay
DefaultAssay(K4.single_cell)  <- assay

min_reads = 5
peaks.common <- table(c(rownames(K4.spatial),rownames(K4.single_cell))) == 2
peaks.common <- peaks.common[peaks.common]

features.common.table <- table(c(rownames(K4.spatial)[rowSums(K4.spatial[[assay]]@counts) > min_reads],
                                 rownames(K4.single_cell)[rowSums(K4.single_cell[[assay]]@counts) > min_reads]))

peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K4.single_cell,K4.spatial),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,reference = 2
)

integrated <- IntegrateData(
  anchorset = anchors,
  preserve.order = TRUE
)

integrated <- RunSVD(
  object = integrated,
  n = 50,
  reduction.name = 'integratedLSI'
)

integrated <- RunUMAP(
  object = integrated,
  dims = 2:40,
  reduction = 'integratedLSI')

p1 <- DimPlot(integrated[,integrated$orig.ident != 'spatial'],pt.size=0.1,label=TRUE)
p2 <- DimPlot(integrated[,integrated$orig.ident == 'spatial'],pt.size=0.1,label=TRUE,group.by='YD.clusters')
p1+p2



### Integrate with scCUT&Tag reference data (H3K27me3)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163532
bin=5000
assay = paste0('bins_',bin)

K27.spatial <- readRDS(file='/data/proj/GCB_MB/spatial_cut-tag/results/CT/H3K27me3/clustering/H3K27me3_clustering.Rds')

p1 <- DimPlot(K27.spatial)
p2 <- DimPlot(K27.spatial,group.by = 'YD.clusters')
p1+p2

K27.single_cell <- readRDS(file='/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/H3K27me3/clustering/01.clustering.Rds')

DefaultAssay(K27.spatial)      <- assay
DefaultAssay(K27.single_cell)  <- assay

min_reads = 5
features.common.table <- table(c(rownames(K27.spatial)[rowSums(K27.spatial[[assay]]@counts) > min_reads],
                                 rownames(K27.single_cell)[rowSums(K27.single_cell[[assay]]@counts) > min_reads]))

peaks.use <- names(features.common.table[features.common.table == 2])

anchors <- FindIntegrationAnchors(
  object.list = list(K27.single_cell,K27.spatial),
  anchor.features = peaks.use,
  assay = rep(assay,2),
  k.filter = NA,reference = 1
)

integrated <- IntegrateData(
  anchorset = anchors,
  preserve.order = TRUE
)

integrated <- RunSVD(
  object = integrated,
  n = 50,
  reduction.name = 'integratedLSI'
)

integrated <- RunUMAP(
  object = integrated,
  dims = 2:40,
  reduction = 'integratedLSI')

p1 <- DimPlot(integrated[,integrated$orig.ident != 'spatial'],pt.size=0.1,label=TRUE)
p2 <- DimPlot(integrated[,integrated$orig.ident == 'spatial'],pt.size=0.1,label=TRUE,group.by='YD.clusters')
p1+p2

SpatialDimPlot(K27.spatial,group.by = 'YD.clusters',pt.size.factor = 6)





