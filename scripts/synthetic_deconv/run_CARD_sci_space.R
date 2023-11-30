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
  
  return(list("card" = props_0, "card_sp" = props_sp, 
              "optimal_phi" = OptimalPhi, "ResList" = ResList))
}


###### run CARD deconvolution on sci-space data
deconv_all_slides_card <- function(data_dir, num_slides = 14){
  set.seed(20220801)
  # load ref expression profiles
  c_x <- read.csv(paste0(data_dir, "/ref_avg_norm_count_markers.csv"), row.names = 1) %>%
    as.matrix() %>%
    as(., "sparseMatrix")
  
  res_l <- c()
  runtime <- c()
  for (s in 1:num_slides){
    cat(paste0("Deconvoluting slide ", s, "\n"))
    
    # load spatial count data and coordinates
    c_y <- Matrix::readMM(
      paste0(data_dir, "/slide_", s, "/syn_sp_count_markers.mtx")
    ) %>%
      t()
    meta <- read.csv(
      paste0(data_dir, "/slide_", s, "/meta_grids.csv"),
      row.names = 1
    )
    coords <- meta[,c("Row", "Col")] %>% `colnames<-`(c("x", "y"))
    
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

save_deconv_results_s <- function(res_dir, res_list){
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
  
  for (s in 1:length(res)){
    # save deconv results
    write.table(res[[s]]$card,
                file = paste0(res_dir, "/card_s", s, ".txt"),
                quote = F, row.names = F, col.names = F
    )
    write.table(res[[s]]$card_sp,
                file = paste0(res_dir, "/card_sp_s", s, ".txt"),
                quote = F, row.names = F, col.names = F
    )
  }
}

data_dir <- "/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/sci_space/deconv_inputs/"
res_dir <- "/Users/jysumac/Projects/Smoother_paper/results/synthetic_deconv/sci_space/card/"

res_s_list <- deconv_all_slides_card(data_dir, num_slides = 14)
save_deconv_results_s(res_dir, res_s_list)