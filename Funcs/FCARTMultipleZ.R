## implement FCART function with multiple covariates
fcart_multi <- function(dat1.list, dat2.list, Z_mat, Lt, Ly, types, vs, optns, bws = c(NULL, NULL)) {
  #n <- length(dat.list)
  dx <- length(vs)
  W_s.list <- list()
  Lsi <- lapply(dat1.list, function(x) {x[[paste0(types[1], "DAY")]]/365})
  Lwi <- lapply(dat1.list, '[[', vs[1])
  W_s.list[[1]] <- remove_duplicates(Lsi, Lwi)
  
  Lsi <- lapply(dat2.list, function(x) {x[[paste0(types[2], "DAY")]]/365})
  Lwi <- lapply(dat2.list, '[[', vs[2])
  W_s.list[[2]] <- remove_duplicates(Lsi, Lwi)
  ## univariate FPCA and covariance calculation
  cat("univariate FPCA:", "\n")
  pca_cov <- uFPCACov(W_s.list = W_s.list, optns)
  pca.list <- pca_cov$pca.list
  for (v in 1:dx) {
    pca <- pca.list[[v]]
    pca$phi <- pca$phi[,1:pca$selectK, drop = FALSE]
    pca$lambda <- pca$lambda[1:pca$selectK]
    pca.list[[v]] <- pca
  }
  idx.list <- pca_cov$idx.list
  nt <- length(Lt); ns <- length(pca_cov$W_is)
  nmt <- names(Ly)
  cat("FullPACE scores:", "\n")
  ## FullPACE pc scores
  xifull <- GetFullPCScore(pca.list = pca.list, cr.cov = pca_cov$cr.cov, W_is = pca_cov$W_is, s.list = pca_cov$s.list, idx.list = idx.list)
  impute_res <- list(FCART = NULL, LOCF = NULL)
  cat("functional imputation: ", "\n")
  for (v in 1:dx) {
    Lw <- W_s.list[[v]]$Ly
    Ls <- W_s.list[[v]]$Lt
    nms <- names(Lw)
    pca <- pca.list[[v]]
    pred <- matrix(rep(pca$mu, ns), byrow = TRUE, ncol = length(pca$workGrid)) + xifull[[v]] %*% t(pca$phi)
    #pred <- xifull[[v]] %*% t(pca$phi)
    s.range <- range(pca$workGrid)
    impute_resv = foreach(j = 1:nt, .combine = "rbind") %dopar% {
      tt <- Lt[[j]]
      resj <- rep(NA, length(tt)); resj2 <- resj # copy for LOCF
      if (nmt[j] %in% nms) {
        ss <- Ls[[nmt[j]]]
        resj2 <- sapply(tt, function(y) {last_obs(Lw[[nmt[j]]], ss, y, nearest = FALSE)})
        in.idx <- which(tt >= s.range[1] & tt <= s.range[2])
        na.idx <- intersect(which(is.na(resj)), in.idx)
        subj.idx <- which(nmt[j] == nms[idx.list[[v]]])
        if (length(na.idx) != 0 && length(subj.idx) == 1) {
          resj[na.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[na.idx], mu = pred[subj.idx,])
          #resj2[na.idx] <- sapply(tt[na.idx], function(y) {last_obs(Lw[[nmt[j]]], ss, y, nearest = FALSE)})
        }
      }
      cbind(resj, resj2)
    }
    impute_res[["FCART"]] <- cbind(impute_res[["FCART"]], impute_resv[,1])
    impute_res[["LOCF"]] <- cbind(impute_res[["LOCF"]], impute_resv[,2])
  }
  colnames(impute_res[["FCART"]]) <- vs
  colnames(impute_res[["LOCF"]]) <- vs
  coef_all <- NULL
  for (m in 1:2) {
    imputem <- impute_res[[m]]
    nna.idx <- which(apply(imputem, 1, function(x) {all(!is.na(x))}))
    #x_rescale <- imputem[nna.idx]/scale_factor
    df <- data.frame(y = unlist(Ly)[nna.idx], imputem[nna.idx,], Z_mat[nna.idx,])
    tt <- unlist(Lt)[nna.idx]
    ett <- seq(min(tt), max(tt), length.out = 100)
    if (m==2 & is.null(bws[m])) {
      tvmodel <- tvLM(as.formula(paste0("y~0+", paste(names(df)[-1], collapse = "+"))), data = df, z = tt, ez = ett, cv.block = length(Lt)/10, 
                      tkernel = "Gaussian", est="ll", bw = bb)
    } else {
      tvmodel <- tvLM(as.formula(paste0("y~0+", paste(names(df)[-1], collapse = "+"))), data = df, z = tt, ez = ett, cv.block = length(Lt)/10, 
                      tkernel = "Gaussian", est="ll", bw = bws[m])
    }
    if (m==1) {bb <- tvmodel$bw}
    coef = as.data.frame(tvmodel$coefficients)
    #names(coef)[1] <- "intercept"
    coef$time <- ett; coef$method <- ifelse(m==1, "FCART", "LOCF")
    XX <- as.matrix(df[,-1]); Y <- df$y
    beta <- solve(t(XX) %*% XX) %*% t(XX) %*% Y
    coef$ti_b0 <- beta[1]; coef$ti_b1 <- beta[2]; 
    coef$bw <- tvmodel$bw
    coef_all <- rbind(coef_all, coef)
  }
  return(list(coef_all, pca_cov))
}

remove_duplicates <- function(oldt, oldy) {
  for (i in 1:length(oldy)) {
    Lti <- oldt[[i]]; Lyi <- oldy[[i]]
    if (length(unique(Lti)) != length(Lti)) {
      Lyi <- sapply(unique(Lti), function(x) {mean(Lyi[which(x == Lti)])})
      oldy[[i]] <- Lyi
      oldt[[i]] <- unique(Lti)
    }
  }
  return(list(Ly = oldy, Lt = oldt))
} 


