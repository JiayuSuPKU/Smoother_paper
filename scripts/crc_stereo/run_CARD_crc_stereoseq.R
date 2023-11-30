library(CARD)
library(tidyverse)

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

save_deconv_results <- function(res_dir, res, runtime = NULL) {
  # check result directory
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = T)
  }
  
  # write down best phi value
  file <- file(paste0(res_dir, "/card_logs.txt"), "w")
  write(paste("Optimal phi: ", as.character(res$optimal_phi), "\n"), file = file)

  if (! is.null(runtime)){
    # save runtime information
    elapsed_time <- runtime["elapsed"] %>% as.character()
    user_time <- runtime["user.self"] %>% as.character()
    sys_time <- runtime["sys.self"] %>% as.character()
    
    cat("Elapsed time:", elapsed_time, "seconds\n")

    write(paste("Elapsed time: ", elapsed_time, " seconds"), file = file)
    write(paste("User CPU time: ", user_time, " seconds"), file = file)
    write(paste("System CPU time: ", sys_time, " seconds"), file = file)
  }
  close(file)

  # save deconv results
  write.table(res$card,
              file = paste0(res_dir, "/card_sp0.txt"),
              quote = F, row.names = F, col.names = T
  )
  write.table(res$card_sp,
              file = paste0(res_dir, "/card_sp1.txt"),
              quote = F, row.names = F, col.names = T
  )
}


library(Seurat)
library(SeuratDisk)
data_dir <- "~/Projects/Smoother_paper/data/crc_stereo/"
res_dir <- "~/Projects/Smoother_paper/results/crc_stereo/"

# load sc ref and synthetic spatial data
sc_data <- LoadH5Seurat(paste0(data_dir, "CRC_SingleCell.h5Seurat"))
sc_labels <- read.csv(paste0(data_dir, "scref/cell_type_sub.csv"), row.names = 1)
sc_data <- sc_data[, rownames(sc_labels)]
sc_data <- AddMetaData(sc_data, sc_labels, "cell_type_sub")
sc_data <- AddMetaData(sc_data, "sample1", "sampleInfo")
ct_names <- levels(as.factor(sc_data$cell_type_sub))

sp_data <- LoadH5Seurat(paste0(data_dir, "P19_T.h5Seurat"))
coords <- sp_data@meta.data[,c('x', 'y')]

# joint QC (all experiments)
card_obj <- createCARDObject(
  sc_count = GetAssayData(sc_data, 'counts'),
  sc_meta = sc_data@meta.data,
  spatial_count = GetAssayData(sp_data, 'counts'),
  spatial_location = coords,
  ct.varname = "cell_type_sub",
  ct.select = ct_names,
  sample.varname = "sampleInfo",
  minCountGene = 0,
  minCountSpot = 0
)
  
# extract deconv inputs
inputs <- extract_card_inputs(card_obj)

c_x <- inputs$ref[, ct_names]
lc_x <- log1p(c_x)
print(dim(c_x))
write.csv(c_x, file = paste0(data_dir, "/scref/ref_avg_count_nmc.csv"))
write.csv(lc_x, file = paste0(data_dir, "/scref/ref_avg_logcount_nmc.csv"))

# release memory
rm(sc_data, sp_data)
gc()

# run deconvolution using CARD default reference genes
runtime <- system.time(
  res <- CARD_deconv_from_ref(
    ref_exp = c_x, spatial_count = inputs$sp, spatial_location = coords
  )
)
save_deconv_results(paste0(res_dir, "p19_t/nmc/"), res, runtime)


# run deconvolution using top50 marker genes
c_x <- read.csv(paste0(data_dir, "/scref/ref_avg_count_nm50.csv"), 
                row.names = 1, check.names = F) %>%
  as.matrix()
sp_data <- LoadH5Seurat(paste0(data_dir, "P19_T.h5Seurat"))
c_y <- GetAssayData(sp_data, 'counts')

gene_intersect = intersect(rownames(c_x), rownames(c_y))
c_x <- c_x[gene_intersect, ]
c_y <- c_y[gene_intersect, ]

rm(sp_data)
gc()

# run deconvolution using truncated marker genes
runtime2 <- system.time(
  res2 <- CARD_deconv_from_ref(
    ref_exp = c_x, spatial_count = c_y, spatial_location = coords
  )
)
save_deconv_results(paste0(res_dir, "p19_t/nm50/"), res2, runtime2)
