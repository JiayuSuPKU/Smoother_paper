library(Seurat)
library(tidyverse)

data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/sci_space/"
res_dir <- "/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/sci_space/deconv_inputs/"


# load sci-space raw count data
counts <- Matrix::readMM(
  paste0(data_dir, "GSE166692_sciSpace_count_matrix.mtx.gz")
)

# load annotations
meta_cells <- read.csv(
  paste0(data_dir, "GSE166692_sciSpace_cell_metadata.tsv.gz"),
  sep = "\t", row.names = 1
)
meta_genes <- read.csv(
  paste0(data_dir, "GSE166692_sciSpace_gene_metadata.tsv.gz"),
  sep = "\t", row.names = 1
)

colnames(counts) <- meta_cells$Cell
rownames(counts) <- meta_genes$id

# preprocess data
sci_space <- CreateSeuratObject(
  counts = counts, meta.data = meta_cells,
  min.cells = 0, min.features = 0
)
sci_space <- NormalizeData(sci_space,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)
sci_space@meta.data$cell_type = sci_space$final_cluster_label %>%
  gsub(" ", "", .)

# select celltypes that are relatively abundant in all slides (>= 100 cells)
ct_count_slide <- table(sci_space$cell_type, sci_space$slide_id)
keep_cts <- (((ct_count_slide >= 100) %>% rowSums()) == ncol(ct_count_slide)) %>%
  which() %>%
  names()

# load marker gene list provided in the original publication (data file S1)
# https://www.science.org/doi/10.1126/science.abb9536
markers <- read_tsv(paste0(data_dir, "sci_space_markers.txt")) %>%
  mutate(cell_type = gsub(' ', '', manual_annotation)) %>%
  filter(cell_type %in% keep_cts)

selected_markers <- unique(
  markers$gene[markers$gene %in% meta_genes$gene_short_name]
)
selected_gids <- meta_genes$id[
  sapply(selected_markers, function(x) {
	which(meta_genes$gene_short_name == x)
  })
]

# generate reference expression files
# raw counts
x_raw_c <- AverageExpression(
  sci_space,
  features = selected_gids,
  group.by = "cell_type", slot = "count"
)$RNA
x_raw_c_ct <- x_raw_c[, keep_cts]
write.csv(
  x_raw_c,
  paste0(res_dir, "/ref_avg_raw_count_markers_all_cts.csv"),
  quote = F
)
write.csv(
  x_raw_c_ct,
  paste0(res_dir, "/ref_avg_raw_count_markers.csv"),
  quote = F
)

# depth-normalized counts
x_norm_c <- AverageExpression(
  sci_space,
  features = selected_gids,
  group.by = "cell_type", slot = "data"
)$RNA
x_norm_c_ct <- x_norm_c[, keep_cts]
write.csv(
  x_norm_c,
  paste0(res_dir, "/ref_avg_norm_count_markers_all_cts.csv"),
  quote = F
)
write.csv(
  x_norm_c_ct,
  paste0(res_dir, "/ref_avg_norm_count_markers.csv"),
  quote = F
)


# assign cells to spots
cell2spot <- function(meta_cells, k = 4,
					  shared_prop = 0.5, mu_n_cells = 10) {
  set.seed(20220801)

  # extract spot info from meta data
  pos_spots <- meta_cells %>%
	group_by(Row, Col) %>%
	summarise(n_cells = n(), .groups = "drop") %>%
	mutate(id = paste(Row, Col, sep = "_")) %>%
	column_to_rownames(., var = "id")

  # get raw cell-spot assignment
  cell2spot <- apply(pos_spots, MARGIN = 1, function(spot) {
	(meta_cells$Row == spot["Row"]) & (meta_cells$Col == spot["Col"])
  }) * 1.0
  rownames(cell2spot) <- rownames(meta_cells)
  colnames(cell2spot) <- rownames(pos_spots)

  # calculate spot-wise distance
  pos <- as.matrix(pos_spots[, c("Row", "Col")])
  dist_sq <- outer(rowSums(pos^2), rowSums(pos^2), "+") -
	2 * tcrossprod(pos, pos)
  dist_sq[dist_sq < 0] <- 0 # for numeric stability

  # calculate binary weight matrix of k-nearest neighbors
  weights <- apply(dist_sq, 1, function(x) {
	thresholds <- unique(sort(x, decreasing = F))[k]
	(x <= thresholds) * 1.0 %>% return()
  })
  # make it symmetric (mutual neighbors)
  weights <- weights * t(weights) * shared_prop
  diag(weights) <- 1.0

  # borrow cells from neighbors
  cell2spot_sm <- cell2spot %*% weights

  # re-sample cells per spot
  cell2spot_rsp <- apply(cell2spot_sm, 2, function(col) { # iterate over spots
	# initialize with no cell selected
	col_sp <- rep(0, length(col))

	# sample cells and control for the expected total number of cells per spot
	n_cells_col <- rpois(1, rgamma(1, mu_n_cells)) # number of cells at the spot
	pool <- which(col > 0) # cells to choose from
	if (length(pool) == 1) { # only one cell at the spot
	  col_sp[pool] <- n_cells_col
	  return(col_sp)
	}
	# sampling with replacement
	cells <- sample(which(col > 0), n_cells_col,
	  replace = T,
	  prob = col[which(col > 0)]
	)
	for (c in cells) { # count occurrence
	  col_sp[c] <- col_sp[c] + 1
	}
	return(col_sp)
  })
  rownames(cell2spot_rsp) <- rownames(meta_cells)
  colnames(cell2spot_rsp) <- rownames(pos_spots)

  # one-hot encode cell type label
  cell_oh_anno <- meta_cells %>%
	select(c(Cell, final_cluster_label)) %>%
	mutate(value = 1) %>%
	pivot_wider(
	  names_from = final_cluster_label, values_from = value, values_fill = 0
	) %>%
	column_to_rownames(var = "Cell") %>%
	as.matrix()
  colnames(cell_oh_anno) <- gsub(" ", "", colnames(cell_oh_anno))

  # generate spot-level annotations before smoothing and resampling
  meta_grids <- meta_cells %>%
	mutate(ct = gsub(" ", "", final_cluster_label)) %>%
	group_by(Row, Col, ct) %>%
	count() %>%
	spread(ct, n, 0) %>%
	ungroup() %>%
	mutate(total = rowSums(select(., -c(Row, Col)))) %>%
	mutate(id = paste(Row, Col, sep = "_")) %>%
	column_to_rownames(., var = "id")
  meta_grids <- meta_grids[rownames(pos_spots), ]

  # spot-level annotations after smoothing and resampling
  meta_grids_rsp <- t(cell2spot_rsp[rownames(cell_oh_anno), ]) %*% cell_oh_anno %>%
	data.frame() %>%
	mutate(
	  total = rowSums(.),
	  Row = pos_spots$Row, Col = pos_spots$Col
	)
  rownames(meta_grids_rsp) <- rownames(pos_spots)

  return(list(
	"cell2spot_raw" = cell2spot, "cell2spot_resampled" = cell2spot_rsp,
	"meta_grids_raw" = meta_grids, "meta_grids_resampled" = meta_grids_rsp
  ))
}

generate_st <- function(slide_id, data = sci_space,
						celltypes = keep_cts, genes = selected_gids,
						resample = TRUE, k = 4, shared_prop = 0.5,
						mu_n_cells = 10, max_umi = 5000, return_data = FALSE,
						res_dir = NULL) {
  set.seed(20220801)

  # select cells in the chosen slide and are of selected cell types
  slide_ct_ind <- (data$cell_type %in% celltypes) &
	(data$slide_id == paste0("Slide ", slide_id))
  meta_cells_s <- data@meta.data[slide_ct_ind, ]

  # assign cells to spots
  res <- cell2spot(meta_cells_s,
	k = k, shared_prop = shared_prop,
	mu_n_cells = mu_n_cells
  )

  if (resample) { # use balanced assignment generated by resampling
	meta_grids <- res$meta_grids_resampled
	cell2spot <- res$cell2spot_resampled
  } else { # use naive pooling of raw data
	meta_grids <- res$meta_grids_raw
	cell2spot <- res$cell2spot_raw
  }

  # generate st counts
  st_counts <- GetAssayData(data[, slide_ct_ind], slot = "count") %*% cell2spot
  colnames(st_counts) <- rownames(meta_grids)

  # downsampling
  st_counts <- SampleUMI(st_counts, max.umi = max_umi, upsample = FALSE, verbose = FALSE)

  # add gamma noise (hyperparameter from the cell2location simulation)
  gene_gamma_noise <- rgamma(shape = 0.5, rate = 0.7, n = nrow(st_counts)) / 2
  st_mean_counts_rg <- st_counts * gene_gamma_noise
  st_counts_rg <- apply(st_mean_counts_rg, MARGIN = 2, function(x) {
	rpois(length(x), lambda = x)
  })
  rownames(st_counts_rg) <- rownames(st_counts)
  colnames(st_counts_rg) <- colnames(st_counts)

  genes <- intersect(genes, rownames(st_counts))
  # n_spots x n_genes
  sp_counts <- st_counts_rg[genes, ] %>% t() %>% as(., "sparseMatrix")

  if (return_data) {
	return(list(
	  "meta_grids" = meta_grids, "cell2spot" = cell2spot,
	  "sp_counts" = sp_counts
	))
  } else {
	stopifnot(!is.null(res_dir))

	# save data
	out_dir <- paste0(res_dir, "/slide_", slide_id)
	if (!dir.exists(out_dir)) {
	  dir.create(out_dir, recursive = T)
	}

	write.csv(meta_grids, paste0(out_dir, "/meta_grids.csv"), quote = F)
	Matrix::writeMM(sp_counts, paste0(out_dir, "/syn_sp_count_markers.mtx"))
  }
}

# generate data for all slides
for (i in 1:14) {
  cat(paste0("Generate spatial data for slide ", i, "\n"))
  generate_st(slide_id = i, data = sci_space, return_data = F, res_dir = res_dir)
}
