library(CARD)
library(tidyverse)

###### prepare CARD inputs for benchmark (donut data)
library(reticulate)
use_condaenv(condaenv = "smoother", conda = "/Users/jysumac/miniforge3/etc/conda")
sc <- import("scanpy")

load_data <- function(data_dir){
  # load single cell reference data
  sc_ref <- sc$read(paste0(data_dir,"/paired_sc_adata.h5ad"))
  sc_meta <- sc_ref$obs["annotation_1"] %>% `colnames<-`("cellType")
  sc_meta$sampleInfo <- "sample1"
  sc_counts <- sc_ref$layers[["counts"]] %>%
    BiocGenerics::t() %>%
    `colnames<-`(., sc_ref$obs_names$values %>% as.character()) %>%
    `rownames<-`(., sc_ref$var_names$values %>% as.character())
  rm(sc_ref)
  gc()
  
  # load spatial data
  sp_syn <- sc$read(paste0(data_dir, "/synthetic_sp_adata.h5ad"))
  sp_counts <- sp_syn$X %>%
    BiocGenerics::t() %>%
    `colnames<-`(., sp_syn$obs_names$values %>% as.character()) %>%
    `rownames<-`(., sp_syn$var_names$values %>% as.character())
  coords <- sp_syn$obsm["X_spatial"][, , 1] %>%
    data.frame() %>%
    `colnames<-`(., c("x", "y")) %>%
    `rownames<-`(., sp_syn$obs_names$values %>% as.character())
  rm(sp_syn)
  gc()
  
  return(list('sc_counts' = sc_counts, 'sc_meta' = sc_meta,
              'sp_counts' = sp_counts, 'coords' = coords))
}

extract_card_inputs <- function(card_obj){
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
  Xinput = spatial_count
  rm(spatial_count)
  B = Basis
  rm(Basis)
  ##### match the common gene names
  Xinput = Xinput[order(rownames(Xinput)),]
  B = B[order(rownames(B)),]
  B = B[rownames(B) %in% common,]
  Xinput = Xinput[rownames(Xinput) %in% common,]
  ##### filter out non expressed genes again
  Xinput = Xinput[rowSums(Xinput) > 0,]
  
  return(list('ref'=B, 'sp'=Xinput))
}

calc_and_save_card_inputs <- function(adata_dir, card_data_dir, num_exp = 10) {
  # check result directory
  out_dir = paste0(card_data_dir, "/deconv_inputs/")
  if (! dir.exists(out_dir)){
    dir.create(out_dir, recursive = T)
  }
  
  # load sc ref and synthetic spatial data
  cat("## Loading H5 data! ...\n")
  data <- load_data(adata_dir)

  # load cell type names
  ab_true <- read.csv(paste0(adata_dir, "/celltype_abundances.csv"), row.names = 1) 
  cell_types <- colnames(ab_true)
  
  # joint QC (all experiments)
  card_obj <- createCARDObject(
    sc_count = data$sc_counts,
    sc_meta = data$sc_meta,
    spatial_count = data$sp_counts,
    spatial_location = data$coords,
    ct.varname = "cellType",
    ct.select = unique(data$sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 0,
    minCountSpot = 0
  )
  
  # extract deconv inputs
  inputs <- extract_card_inputs(card_obj)
  
  cat("## Saving input data! ...\n")
  # copy the cell type abundances to the new directory
  write.csv(ab_true, file = paste0(card_data_dir, "/celltype_abundances.csv"), quote = F)
  # copy the cell type assignment matrix (with high low density anno) to the new directory
  file.copy(paste0(adata_dir, "/celltype_zone_assignment.csv"), paste0(card_data_dir, "/celltype_zone_assignment.csv"))

  c_x <- inputs$ref[, cell_types]
  write.csv(c_x, file = paste0(out_dir, "/ref_avg_norm_count_markers.csv"))
  
  for (exp in 0:(num_exp - 1)){
    exp_inds <- (2500 * exp + 1):(2500 + 2500 * exp)
    c_y <- inputs$sp[, exp_inds] %>% t()
    Matrix::writeMM(c_y, paste0(out_dir, "/syn_sp_count_markers_exp", exp, ".mtx"))
    coords_exp <- data$coords[exp_inds, ]
    write.csv(coords_exp, file = paste0(
      out_dir, "/coords_exp", exp, ".csv"))
  }
}

###### prepare CARD inputs for benchmark (donut data)
in_dir <- "/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/donut/"
cat("## Extracting CARD-selected reference marker genes ...\n")
# extract and save inputs as csv
for (sm in c('0', '0.1', '0.5')){
  data_dir <- paste0(
    in_dir,
    "ne10_rz15_nm0_sm",
    sm, '/'
  )
  card_data_dir <- paste0(
    in_dir,
    "ne10_rz15_nmc_sm",
    sm, '/'
  )
  calc_and_save_card_inputs(data_dir, card_data_dir)
}

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
  
  return(list("card" = props_0, "card_sp" = props_sp, 
              "optimal_phi" = OptimalPhi, "ResList" = ResList))
}

deconv_all_exp_card <- function(data_dir, num_exp = 10){
  set.seed(20220801)
  # load ref expression profiles
  ab_true <- read.csv(paste0(data_dir, "/celltype_abundances.csv"), row.names = 1)
  c_x <- read.csv(paste0(data_dir, "/deconv_inputs/ref_avg_norm_count_markers.csv"), 
                  row.names = 1) %>%
    as.matrix() %>%
    as(., "sparseMatrix")
  # make sure cell types are in the correct order
  c_x <- c_x[, colnames(ab_true)]
  
  res_l <- c()
  runtime <- c()
  for (exp in 0:(num_exp - 1)){
    cat(paste0("Deconvoluting exp ", exp, "\n"))
    
    # load spatial count data and coordinates
    c_y <- Matrix::readMM(
      paste0(data_dir, "/deconv_inputs/syn_sp_count_markers_exp", exp, ".mtx")
    ) %>%
      t()
    coords <- read.csv(
      paste0(data_dir, "/deconv_inputs/coords_exp", exp, ".csv"),
      row.names = 1
    )
    
    t_start <- Sys.time()
    # run deconvolution
    res <- CARD_deconv_from_ref(
      ref_exp = c_x, spatial_count = c_y, spatial_location = coords
    )
    t_end <- Sys.time()
    
    # save results
    res_l <- c(res_l, list(res))
    runtime <- c(runtime, difftime(t_end, t_start, units = 'secs'))
    
    # clear memory
    rm(c_y, coords, res)
    gc()
  }
  
  return(list('result'=res_l, 'runtime'=runtime))
}

save_deconv_results <- function(res_dir, res_list){
  # check result directory
  if (! dir.exists(res_dir)){
    dir.create(res_dir, recursive = T)
  }
  
  # save runtimes
  runtime <- res_list$runtime %>% data.frame() %>%
    `colnames<-`('Runtime (secs)')
  write.csv(runtime, file = paste0(res_dir, "/runtime.csv"), quote = F)
  
  res <- res_list$result
  # write down best phi value
  optim_phi = sapply(res, function(x) x$optimal_phi) %>% 
    data.frame() %>% `colnames<-`('Optimal phi')
  write.csv(optim_phi, file = paste0(res_dir, "/optim_phi.csv"), quote = F)
  
  for (exp in 1:length(res)){
    # save deconv results
    write.table(res[[exp]]$card,
                file = paste0(res_dir, "/card_e", exp - 1, ".txt"),
                quote = F, row.names = F, col.names = F
    )
    write.table(res[[exp]]$card_sp,
                file = paste0(res_dir, "/card_sp_e", exp - 1, ".txt"),
                quote = F, row.names = F, col.names = F
    )
  }
}

cat("## Run CARD deconvolution ...\n")
in_dir <- "/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/donut/"
out_dir <- "/Users/jysumac/Projects/Smoother_paper/results/synthetic_deconv/donut/"

# deconvolute all data
for (nm in c('20', '50', '0', 'c')){
  for (sm in c('0', '0.1', '0.5')){
    data_dir <- paste0(
      in_dir,
      sprintf(
        "ne10_rz15_nm%s_sm%s/",
        nm, sm
      )
    )
    res_dir <- paste0(
      out_dir,
      sprintf(
        "ne10_rz15_nm%s_sm%s/card/",
        nm, sm
      )
    )
    print(data_dir)
    
    res_exp_list <- deconv_all_exp_card(data_dir, num_exp = 10)
    save_deconv_results(res_dir, res_exp_list)
  }
}