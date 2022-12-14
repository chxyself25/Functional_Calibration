## asynchronous longitudinal regression for one X case
## wrap up FPCA, calibration, and regression processses
## input: two functional observations including both Y and t in a list format
FMlassoX1 <- function(Y, t, Lw, Ls, optns = NULL) {
  if (is.null(optns)) {
    optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
  }
  n <- length(Lw) # Y and W should have the same number of subjects
  pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
  pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi[, 1:pca$selectK, drop = FALSE])
  s.range <- range(pca$workGrid)
  #impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
  impute_res = unlist(lapply(1:n, function(j) {
    tt <- t[j,]
    in.idx <- tt >= s.range[1] & tt <= s.range[2]
    resj <- rep(NA, length(tt))
    resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
    resj
  }))
  nna.idx <- which(!is.na(impute_res))
  X.fm <- cbind(1, impute_res[nna.idx])
  Y.nna <- Y[nna.idx]
  b.fm <- solve(t(X.fm) %*% X.fm) %*% t(X.fm) %*% Y.nna
  # test on beta_x
  sigma2 <- (t(Y.nna)%*%Y.nna - t(Y.nna)%*% X.fm %*% b.fm)/(length(Y.nna)-2)
  #cat(t(X.fm)%*%X.fm, "\n")
  beta.sigma2 <- c(sigma2) *  solve(t(X.fm)%*%X.fm)
  test.stat <- b.fm[2]/sqrt(beta.sigma2[2,2])
  prob <- 2*pt(abs(test.stat), df = length(Y.nna)-2, lower.tail = FALSE)
  return(list(beta = b.fm, beta.sigma2 = diag(beta.sigma2),tval = test.stat, pval = prob,
              pca.sigma2 = pca$sigma2, pca.rho = pca$rho, pca.K = pca$selectK, pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
}


# function that can return the nearest last observation indexes of s
# if there is no last observation, then return the nearest observation
last_obs <- function(ws=NULL, s, tt, nearest = FALSE) {
  idx1 <- which(s <= tt)
  if (length(idx1) == 0) {
    if (!nearest) {return(NA)}
    else {return(ws[which.min(abs(s-tt))])}
  } else {
    return(ws[max(idx1)])
    #return(ws[idx1[which.max(s[idx1])]])
  }
}

## function for last observation carried forward
LOCFlassoX <- function(Y, t, Ws, s) {
  n <- nrow(t)
  m <- ncol(t)
  dx <- ncol(Ws)
  nXhat <- matrix(0, ncol = dx, nrow = m*n)
  for (v in 1:dx) {
    Wsv <- matrix(Ws[,v], ncol = m, nrow = n, byrow = TRUE)
    sv <- s[[v]]
    nXhat[,v] <- c(sapply(1:n, function(i) {
      sapply(t[i,], function(y) {last_obs(Wsv[i,], sv[i,], y, nearest = FALSE)})
      }))
  }
  nna.idx <- which(apply(nXhat, 1, function(x) {all(!is.na(x))}))
  Xmat <- cbind(1, nXhat[nna.idx,])
  Y.nna <- Y[nna.idx]
  beta <- solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% Y.nna
  # test on beta_x
  sigma2 <- (t(Y.nna)%*%Y.nna - t(Y.nna) %*% Xmat %*% beta)/(length(Y.nna)-2)
  beta.sigma2 <- c(sigma2) * solve(t(Xmat)%*%Xmat)
  test.stat <- beta[2]/sqrt(beta.sigma2[2,2])
  prob <- 2*pt(abs(test.stat), df = length(Y.nna)-2, lower.tail = FALSE)
  return(list(beta = beta, beta.sigma2 = diag(beta.sigma2), tval = test.stat, pval = prob))
}

# time-dependent coefficient: only one X
tvFMlassoX1 <- function(Y, t, Lw, Ls, optns = list(), convert = TRUE, bw = NULL) {
  if (is.null(optns)) {
    optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
  }
  n <- length(Lw) # Y and W should have the same number of subjects
  pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
  pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi[, 1:pca$selectK, drop = FALSE])
  s.range <- range(pca$workGrid)
  impute_res = foreach(j = 1:n, .combine = "c") %dopar% {
    tt <- t[j,]
    in.idx <- tt >= s.range[1] & tt <= s.range[2]
    resj <- rep(NA, length(tt))
    resj[in.idx] <- ConvertSupport(pca$workGrid, toGrid = tt[in.idx], mu = pred[j,])
    resj
  }
  nna.idx <- which(!is.na(impute_res))
  tt <- c(t(t))[nna.idx]
  df <- data.frame(y = Y[nna.idx][order(tt)], x = impute_res[nna.idx][order(tt)])
  #tvmodel <- tvLM(y~x, data = df, cv.block = n/10, est="ll")
  if (convert) {
    ttsort <- seq(min(tt), max(tt), length.out = optns$nRegGrid)
    #coef <- ConvertSupport(fromGrid = tt[order(tt)], toGrid = ttsort, phi = tvmodel$coefficients)
  } else {
    ttsort <- tt[order(tt)]
    #coef <- tvmodel$coefficients
  }
  tvmodel <- tvLM(y~x, z = tt[order(tt)], ez = ttsort, data = df, cv.block = n/10, est="ll", bw = bw)
  coef = tvmodel$coefficients
  diff.mat <- cbind(beta_fun(ttsort, intercept = TRUE), beta_fun(ttsort, intercept = FALSE)) - coef
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint, coef= coef, time = ttsort,
              bw = tvmodel$bw, pca.K = pca$selectK, pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
}

tvLOCFlassoX1 <- function(Y, t, Lw, Ls, grid.len = 60, bw = NULL) {
  n <- nrow(t)
  m <- ncol(t)
  impute_res <- foreach(j = 1:n, .combine = "c") %dopar% {
    tt <- t[j,]
    sapply(tt, function(x) {last_obs(Lw[[j]], Ls[[j]], x, nearest = FALSE)})
  }
  nna.idx <- which(!is.na(impute_res))
  tt <- c(t(t))[nna.idx]
  df <- data.frame(y = Y[nna.idx][order(tt)], x = impute_res[nna.idx][order(tt)])
  #tvmodel <- tvLM(y~x, data = df, cv.block = n/10, est="ll")
  if (!is.null(grid.len)) {
    ttsort <- seq(min(tt), max(tt), length.out = grid.len)
    #coef <- ConvertSupport(fromGrid = tt[order(tt)], toGrid = ttsort, phi = tvmodel$coefficients)
  } else {
    ttsort <- tt[order(tt)]
    #coef <- tvmodel$coefficients
  }
  tvmodel <- tvLM(y~x, z = tt[order(tt)], ez = ttsort, data = df, cv.block = n/10, est="ll", bw = bw)
  coef = tvmodel$coefficients
  diff.mat <- cbind(beta_fun(ttsort, intercept = TRUE), beta_fun(ttsort, intercept = FALSE)) - coef
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint, 
              bw = tvmodel$bw, coef = coef, time = ttsort))
}

tvCZFlassoX1 <- function(Y, t, Ws, s, grid.len = 60, bw = NULL) {
  n <- nrow(t)
  m <- ncol(t)
  ## rescale s and t
  s <- s/10; t <- t/10
  data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(s)), W = Ws)
  data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(t)), Y = Y)
  #timepts <- sort(c(t))
  timepts <- seq(max(min(t), min(s)), min(max(t), max(s)), length.out = grid.len)
  ti.res <- foreach(x = timepts, .combine = "rbind") %do% {
    timodel <- tryCatch(asynchTD(data.x, data.y, times = x, kType = "epan", lType = "identity",
             bw = bw, verbose = FALSE), error = function(e) {NA})
    if (any(is.na(timodel))) {
      rep(NA, 3)
    } else {
      c(c(timodel$betaHat), ifelse(is.null(bw), timodel$optBW, bw)) 
    }
  }
  diff.mat <- cbind(beta_fun(timepts, intercept = TRUE), beta_fun(timepts, intercept = FALSE)) - ti.res[,1:2]
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {
    nna.idx <- which(!is.na(x)); trapzRcpp(timepts[nna.idx], abs(x[nna.idx]))})
  mseint <- apply(diff.mat, 2, function(x) {
    nna.idx <- which(!is.na(x)); trapzRcpp(timepts[nna.idx], (x[nna.idx])^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint,
              coef = ti.res[,1:2], time = 10*timepts, bw = ti.res[,3]))
}

## convert support for cross covariance
ConvertSuportcr <- function(fromGrid, toGrid, Cov = NULL) {
  for (i in 1:2) {
    fgrid <- fromGrid[,i]; tgrid <- toGrid[,i]
    # In case the range of toGrid is larger than fromGrid due to numeric error
    #buff <- .Machine$double.eps * max(abs(fromGrid)) * 10
    buff <- 10e-06
    if (abs(tgrid[1] - fgrid[1]) < buff)
      tgrid[1] <- fgrid[1]
    if (abs(tgrid[length(tgrid)] - fgrid[length(fgrid)]) < buff)
      tgrid[length(tgrid)] <- fgrid[length(fgrid)]
    if ( ( fgrid[1] - buff  >  tgrid[1]) || 
         ( fgrid[length(fgrid)] + buff < tgrid[length(tgrid)]) ) {
      stop("Insufficient size of 'fromGrid'.")}
    toGrid[,i] <- tgrid
  }
  gd <- expand.grid(X=toGrid[,1], Y=toGrid[,2])
  ret <- matrix(fdapace:::interp2lin(fromGrid[,1], fromGrid[,2], Cov, gd$X, gd$Y), nrow=length(toGrid[,1]))
  return(ret)
}

tvFunctionlassoX1 <- function(Y, t, Lw, Ls, optnsx = list(), optnsy = list(), cr.bw = c(NULL, NULL)) {
  n <- nrow(t)
  m <- ncol(t)
  # estimate mean function of Y
  Ly <- split(Y, as.factor(rep(1:n, each = m)))
  Lt <- lapply(1:n, function(x) {t[x,]})
  muY <- FPCAMu(Ly, Lt, optnsy)
  # estimate covariance and mean of X
  pcaX <- FPCAIC(Lw, Ls, optnsx, spec = list("AIC", FALSE))
  # find a unified grid
  gs <- max(muY$workGrid[1], pcaX$workGrid[1])
  ge <- min(max(muY$workGrid), max(pcaX$workGrid))
  grid <- seq(gs, ge, length.out = max(optnsx$nRegGrid, optnsy$nRegGrid))
  # transform to the grid
  muy <- ConvertSupport(fromGrid = muY$workGrid, toGrid = grid, mu = muY$muWork)
  mux <- ConvertSupport(fromGrid = pcaX$workGrid, toGrid = grid, mu = pcaX$mu)
  CovXX <- ConvertSupport(fromGrid = pcaX$workGrid, toGrid = grid, mu = diag(pcaX$fittedCov))
  # estimate cross-covariance
  cov <- GetCrCovYX(Ly1 = muY$resid, Lt1 = Lt, Ly2 = pcaX$resid, Lt2 = Ls, 
                    Ymu1 = rep(0, length(unique(unlist(Lt)))), 
                    Ymu2 = rep(0, length(unique(unlist(Ls)))), 
                    bw1 = cr.bw[1], bw2 = cr.bw[2], kern = "gauss")
  Cov <- ConvertSuportcr(fromGrid = cov$smoothGrid, toGrid = cbind(grid, grid), Cov = cov$smoothedCC)
  CovYX <- diag(Cov)
  beta1 <- CovYX/CovXX
  beta0 <- muy - beta1*mux
  diff.mat <- cbind(beta_fun(grid, intercept = TRUE), beta_fun(grid, intercept = FALSE)) - cbind(beta0, beta1)
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(grid, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(grid, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint,
              muy.bw = muY$bwMu, pca.mubw = pcaX$bwMu, pca.covbw = pcaX$bwCov, cr.bw = cov$bw,
              coef = cbind(beta0, beta1), time = grid, bw = cov$bw))
}

tvOracle <- function(Y, t, X, grid.len = 60, bw = NULL) {
  tt <- c(t(t))
  df <- data.frame(y = Y[order(tt)], x = X[order(tt)])
  #tvmodel <- tvLM(y~x, data = df, cv.block = n/10, est="ll")
  if (!is.null(grid.len)) {
    ttsort <- seq(min(tt), max(tt), length.out = grid.len)
    #coef <- ConvertSupport(fromGrid = tt[order(tt)], toGrid = ttsort, phi = tvmodel$coefficients)
  } else {
    ttsort <- tt[order(tt)]
    #coef <- tvmodel$coefficients
  }
  tvmodel <- tvLM(y~x, z = tt[order(tt)], ez = ttsort, data = df, cv.block = n/10, est="ll", bw = bw)
  coef = tvmodel$coefficients
  diff.mat <- cbind(beta_fun(ttsort, intercept = TRUE), beta_fun(ttsort, intercept = FALSE)) - coef
  mae <- apply(diff.mat, 2, function(x) {mean(abs(x), na.rm = TRUE)})
  mse <- apply(diff.mat, 2, function(x) {mean(x^2, na.rm = TRUE)})
  maeint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, abs(x))})
  mseint <- apply(diff.mat, 2, function(x) {trapzRcpp(ttsort, x^2)})
  return(list(mae = mae, mse = mse, maeint = maeint, mseint = mseint, coef= coef, time = ttsort,
              bw = tvmodel$bw))
}
############# old functions using Fan's 2 step method ###############
## time-dependent coefficient: only one X 
# tdFMlassoX1 <- function(Ly, t, Lw, Ls, optns = NULL) {
#   if (is.null(optns)) {
#     optns <- list(dataType = "Sparse", nRegGrid = 51, methodBwMu = "GCV", methodBwCov = "GCV")
#   }
#   n <- length(Lw) # Y and W should have the same number of subjects
#   pca <- FPCAIC(Lw, Ls, optns, spec = list("AIC", FALSE))
#   pred <- matrix(rep(pca$mu, n), byrow = TRUE, ncol = length(pca$workGrid)) + pca$xiEst %*% t(pca$phi)
#   s.range <- range(pca$workGrid)
#   tt <- sort(unique(c(t)))
#   in.idx <- tt >= s.range[1] & tt <= s.range[2]
#   tt <- tt[in.idx]; m <- length(tt)
#   impute_res <- lapply(1:n, function(j) {
#     ConvertSupport(pca$workGrid, toGrid = tt, mu = pred[j,])
#   })
#   betat <- foreach(j = 1:m, .combine = "rbind") %dopar% {
#     tj <- tt[j]
#     xj <- sapply(impute_res, '[', j)
#     yj <- sapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         Ly[[x]][which(t[x,]==tj)]
#         } else {
#         NA
#         }
#     })
#     nna.idx <- which(!is.na(yj))
#     yj <- yj[nna.idx]
#     Xj <- cbind(1, xj[nna.idx])
#     bj <- tryCatch(solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj, error = function(e) {rep(NA, 2)})
#     c(bj)
#   }
#   return(list(beta = betat, time = tt, pca.sigma2 = pca$sigma2, pca.rho = pca$rho, pca.K = pca$selectK, 
#               pca.mubw = pca$bwMu, pca.covbw = pca$bwCov))
# }
# 
# ## function for last observation carried forward
# tdLOCFlassoX <- function(Ly, t, Ws, s) {
#   n <- nrow(t)
#   m <- ncol(t)
#   dx <- ncol(Ws)
#   nXhat <- matrix(0, ncol = dx, nrow = m*n)
#   for (v in 1:dx) {
#     Wsv <- matrix(Ws[,v], ncol = m, nrow = n, byrow = TRUE)
#     sv <- s[[v]]
#     nXhat[,v] <- c(sapply(1:n, function(i) {
#       sapply(t[i,], function(y) {last_obs(Wsv[i,], sv[i,], y, nearest = FALSE)})
#     }))
#   }
#   tt <- sort(unique(c(t)))
#   betat <- foreach(j = 1:length(tt), .combine = "rbind") %dopar% {
#     tj = tt[j]
#     yj <- unlist(lapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         Ly[[x]][which(t[x,]==tj)]
#       }
#     }))
#     xj <- lapply(1:n, function(x) {
#       if (tj %in% t[x,]) {
#         nXhat[which(t[x,]==tj)+m*(x-1),,drop=FALSE]
#       }
#     })
#     xj <- do.call("rbind", xj)
#     nna.idx <- which(apply(xj, 1, function(x) {all(!is.na(x))}))
#     Xj <- cbind(1, xj[nna.idx,,drop=FALSE])
#     yj <- yj[nna.idx]
#     #bj <- solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj
#     bj <- tryCatch(solve(t(Xj) %*% Xj) %*% t(Xj) %*% yj, error = function(e) {rep(NA, 2)})
#     c(bj) 
#   }
#   return(list(beta = betat, time = tt))
# }
# 
