library(CARD)
library(tidyverse)

###### prepare CARD inputs for benchmark
library(reticulate)
use_condaenv(condaenv = "smoother", conda = "/Users/jysumac/miniforge3/etc/conda")
sc <- import("scanpy")
ad <- import("anndata")

# load 10x Visium spatial data
load_sp_data <- function(data_dir, sample_id){
  sp_dir <- paste0(data_dir, '/ST/', sample_id, '/')
  h5_dir <- sprintf("%s/%s_filtered_feature_bc_matrix.h5", sp_dir, sample_id)
  sp_data <- sc$read_10x_h5(h5_dir)
  ad$AnnData$var_names_make_unique(sp_data)
  
  sp_counts <- sp_data$X %>%
    BiocGenerics::t() %>%
    `colnames<-`(., sp_data$obs_names$values %>% as.character()) %>%
    `rownames<-`(., sp_data$var_names$values %>% as.character())
  
  coords_dir <- sprintf("%s/spatial/tissue_positions_list.csv", sp_dir)
  coords <- read.csv(coords_dir, header = F, row.names = 1) %>%
    `colnames<-`(., c(
      "in_tissue",
      "array_row",
      "array_col",
      "x",
      "y"
    ))
  coords <- coords[colnames(sp_counts), c('x', 'y')]
  
  rm(sp_data)
  gc()
  
  return(list('sp_counts' = sp_counts, 'coords' = coords))
}

# load sc reference data
load_sc_data <- function(data_dir) {
  # load single cell reference data
  sc_ref <- sc$read(paste0(data_dir, "/scref/scref_with_ct_mks.h5ad"))
  sc_meta <- sc_ref$obs["cell_type"] %>% `colnames<-`("cellType")
  sc_meta$sampleInfo <- "sample1"
  sc_counts <- sc_ref$layers['raw'] %>%
    BiocGenerics::t() %>%
    `colnames<-`(., sc_ref$obs_names$values %>% as.character()) %>%
    `rownames<-`(., sc_ref$var_names$values %>% as.character())
  rm(sc_ref)
  gc()

  return(list("sc_counts" = sc_counts, "sc_meta" = sc_meta))
}

# extract reference and sp count matrix from a card object
extract_card_inputs <- function(card_obj) {
  # extract info
  ct.select <- card_obj@info_parameters$ct.select
  ct.varname <- card_obj@info_parameters$ct.varname
  sample.varname <- card_obj@info_parameters$sample.varname
  # extract data
  sc_eset <- card_obj@sc_eset
  # calculate avg expression profile
  Basis_ref <- createscRef(sc_eset, ct.select, ct.varname, sample.varname)
  Basis <- Basis_ref$basis
  Basis <- Basis[, colnames(Basis) %in% ct.select]
  Basis <- Basis[, match(ct.select, colnames(Basis))]
  spatial_count <- card_obj@spatial_countMat

  # filter genes
  commonGene <- intersect(rownames(spatial_count), rownames(Basis))
  #### remove mitochondrial and ribosomal genes
  commonGene <- commonGene[!(commonGene %in% commonGene[grep("mt-", commonGene)])]
  common <- selectInfo(Basis, sc_eset, commonGene, ct.select, ct.varname)
  Xinput <- spatial_count
  rm(spatial_count)
  B <- Basis
  rm(Basis)
  ##### match the common gene names
  Xinput <- Xinput[order(rownames(Xinput)), ]
  B <- B[order(rownames(B)), ]
  B <- B[rownames(B) %in% common, ]
  Xinput <- Xinput[rownames(Xinput) %in% common, ]
  ##### filter out non expressed genes again
  Xinput <- Xinput[rowSums(Xinput) > 0, ]

  return(list("ref" = B, "sp" = Xinput))
}

# calculate and save inputs used in card
calc_and_save_card_inputs <- function(data_dir, sample_id) {
  # check result directory
  out_dir <- paste0(data_dir, "/scref/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }

  # load sc ref and synthetic spatial data
  cat("## Loading H5 data! ...\n")
  sc_data <- load_sc_data(data_dir)
  sp_data <- load_sp_data(data_dir, sample_id)
  
  ct_names <- read.csv(paste0(data_dir, "/scref/ref_avg_count_nm0.csv"),
    row.names = 1, nrows = 1
  ) %>% colnames()

  # joint QC (all experiments)
  card_obj <- createCARDObject(
    sc_count = sc_data$sc_counts,
    sc_meta = sc_data$sc_meta,
    spatial_count = sp_data$sp_counts,
    spatial_location = sp_data$coords,
    ct.varname = "cellType",
    ct.select = unique(sc_data$sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 0,
    minCountSpot = 0
  )

  # extract deconv inputs
  inputs <- extract_card_inputs(card_obj)

  cat("## Saving input data! ...\n")
  c_x <- inputs$ref[, ct_names]
  lc_x <- log1p(c_x)
  write.csv(c_x, file = paste0(out_dir, "/ref_avg_count_nmc.csv"))
  write.csv(lc_x, file = paste0(out_dir, "/ref_avg_logcount_nmc.csv"))
}

# extract and save inputs as csv
print("## Calculating and saving CARD inputs for mouse brain visium data")
data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/mouse_brain_visium/"
sample_id <- "ST8059048"
calc_and_save_card_inputs(data_dir, sample_id)


###### run CARD deconvolution on donut data
CARD_deconv_from_ref <- function(ref_exp, spatial_count, spatial_location,
                                 isigma = 0.1, epsilon = 1e-04, max_iter = 1000,
                                 init_with_lr_sol = FALSE, inv_cov = NULL) {
  # prepare inputs
  ct.select <- colnames(ref_exp)
  B <- ref_exp
  Xinput <- spatial_count
  colsumvec <- colSums(Xinput) + 1e-10
  Xinput_norm <- sweep(Xinput, 2, colsumvec, "/")

  # scale coordinates
  norm_cords <- spatial_location[, c("x", "y")]
  norm_cords$x <- norm_cords$x - min(norm_cords$x)
  norm_cords$y <- norm_cords$y - min(norm_cords$y)
  scaleFactor <- max(norm_cords$x, norm_cords$y)
  norm_cords$x <- norm_cords$x / scaleFactor
  norm_cords$y <- norm_cords$y / scaleFactor

  ED <- fields::rdist(as.matrix(norm_cords)) ## Euclidean distance matrix

  ##### initialize the proportion matrix
  if (init_with_lr_sol) {
    Vint1 <- (solve(t(B) %*% B) %*% t(B) %*% Xinput_norm) %>% t()
  } else {
    set.seed(20200107)
    Vint1 <- MCMCpack::rdirichlet(ncol(Xinput_norm), rep(10, ncol(B)))
    mean_X <- mean(Xinput_norm)
    mean_B <- mean(B)
    Xinput_norm <- Xinput_norm * 1e-01 / mean_X
    B <- B * 1e-01 / mean_B
  }
  Vint1 <- as.matrix(Vint1)
  colnames(Vint1) <- colnames(B)
  rownames(Vint1) <- colnames(Xinput_norm)

  ###### parameters that need to be set
  isigma <- isigma #### construct Gaussian kernel with the default scale /length parameter to be 0.1
  kernel_mat <- exp(-ED^2 / (2 * isigma^2))
  diag(kernel_mat) <- 0

  epsilon <- epsilon #### convergence epsion
  phi <- c(0, 0.3, 0.6, 0.99) #### grided values for phi

  rm(ED)
  rm(Xinput)
  rm(norm_cords)
  gc()
  ResList <- list()
  Obj <- c()
  for (iphi in 1:length(phi)) {
    res <- CARDref(
      XinputIn = as.matrix(Xinput_norm),
      UIn = as.matrix(B),
      WIn = kernel_mat,
      phiIn = phi[iphi],
      max_iterIn = max_iter,
      epsilonIn = epsilon,
      initV = Vint1,
      initb = rep(0, ncol(B)),
      initSigma_e2 = 0.1,
      initLambda = rep(10, length(ct.select))
    )
    rownames(res$V) <- colnames(Xinput_norm)
    colnames(res$V) <- colnames(B)
    ResList[[iphi]] <- res
    Obj <- c(Obj, res$Obj)
  }
  Optimal <- which(Obj == max(Obj))
  Optimal <- Optimal[length(Optimal)] #### just in case if there are two equal objective function values
  OptimalPhi <- phi[Optimal]
  OptimalRes <- ResList[[Optimal]]

  props_0 <- sweep(ResList[[1]]$V, 1, rowSums(ResList[[1]]$V), "/")
  props_sp <- sweep(OptimalRes$V, 1, rowSums(OptimalRes$V), "/")

  return(list(
    "card" = props_0, "card_sp" = props_sp,
    "optimal_phi" = OptimalPhi, "ResList" = ResList
  ))
}


deconv_card <- function(data_dir, sample_id, nm = '20') {
  # load ref expression profiles
  c_x <- read.csv(sprintf("%s/scref/ref_avg_count_nm%s.csv", data_dir, nm), 
                  header = T, row.names = 1) %>%
    as.matrix() %>%
    as(., "sparseMatrix")
  
  # load spatial counts and coordinates
  sp_data <- load_sp_data(data_dir, sample_id)

  cat(paste0("Deconvoluting sample ", sample_id, "\n"))

  # run deconvolution
  res <- CARD_deconv_from_ref(
    ref_exp = c_x,
    spatial_count = sp_data$sp_counts[rownames(c_x),],
    spatial_location = sp_data$coords
  )

  return(res)
}

save_deconv_results <- function(res_dir, res_list, s_id_list) {
  
  stopifnot(length(res_list) == length(s_id_list))
  
  # check result directory
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = T)
  }
  
  # write down best phi value
  optim_phi <- sapply(res_list, function(x) x$optimal_phi) %>%
    data.frame() %>%
    `colnames<-`("Optim_phi")
  write.csv(optim_phi, file = paste0(res_dir, "/optim_phi.csv"), quote = F)
  
  for (i in 1:length(res_list)) {
    # check result directory
    sample_dir <- paste0(res_dir, "/", s_id_list[i], "/")
    if (!dir.exists(sample_dir)) {
      dir.create(sample_dir, recursive = T)
    }
    
    # save deconv results
    write.table(res_list[[i]]$card,
                file = paste0(sample_dir, "/card_sp0.txt"),
                quote = F, row.names = F, col.names = F
    )
    write.table(res_list[[i]]$card_sp,
                file = paste0(sample_dir, "/card_sp1.txt"),
                quote = F, row.names = F, col.names = F
    )
  }
}

# run CARD deconvolution
print("## Running CARD deconvolution on mouse brain visium data")

# data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/mouse_brain_visium/"
# sample_id <- "ST8059048"
# res <- deconv_card(data_dir, sample_id, nm = 'c')
# res_dir <- "/Users/jysumac/Projects/Smoother_paper/results/mouse_brain_visium/nmc/"
# save_deconv_results(res_dir, list(res), c(sample_id))

data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/mouse_brain_visium/"
sample_id <- "ST8059048"

for (nm in c("20", "50", "0", "c")) {
  res <- deconv_card(data_dir, sample_id, nm = nm)
  res_dir <- paste0(
    "/Users/jysumac/Projects/Smoother_paper/results/mouse_brain_visium/nm", nm, "/"
  )
  save_deconv_results(res_dir, list(res), c(sample_id))
}
