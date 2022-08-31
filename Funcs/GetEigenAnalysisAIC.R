## modified eigenfunction decomposition function
GetEigenAnalysisAIC = function(y, t, smoothCov, workGrid, obsGrid, sigma2, muWork = NULL, muObs, 
                               spec, optns) {
  maxK <- optns$maxK # maximum number of principle components to consider, default: 20
  verbose <- optns$verbose
  n <- length(y); N <- length(unlist(y))
  method <- spec[[1]]; correct <- spec[[2]]
  
  gridSize <- workGrid[2] - workGrid[1]
  numGrids <- nrow(smoothCov)
  
  eig <- eigen(smoothCov)
  
  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]
  
  if (maxK < length(d)) {
    if (optns[['verbose']]) {
      message(sprintf("At most %d number of PC can be selected, thresholded by `maxK` = %d. \n", length(d), maxK)) 
    }
    
    d <- d[1:maxK]
    eigenV <- eigenV[, 1:maxK, drop=FALSE]
  }
  FVE <- cumsum(d)/sum(d) * 100
  
  # normalization
  if (is.null(muWork)) {
    muWork = 1:dim(eigenV)[1]
  }
  phi <- apply(eigenV, 2, function(x) {
    x <- x / sqrt(trapzRcpp(workGrid, x^2)) 
    if ( 0 <= sum(x*muWork) )
      return(x)
    else
      return(-x)
  })
  lambda <- gridSize * d;
  
  # convert phi to obsGrid
  phiObs <- ConvertSupport(workGrid, obsGrid, phi=phi)
  
  # get regularizing scalar for measurement error variance estimate
  if (correct) {
    CovObs <- phiObs %*% diag(lambda) %*% t(phiObs)
    if( length(y) > 2048 ){
      randIndx <- sample( length(y), 2048)
      rho <- GetRho(y[randIndx], t[randIndx], optns, muObs, obsGrid, CovObs, lambda, phiObs, sigma2)
    } else {
      rho <- GetRho(y, t, optns, muObs, obsGrid, CovObs, lambda, phiObs, sigma2)
    }
    if (sigma2 <= rho) {sigma2 <- rho}
  }
  
  ## select number of principle components by two kinds of AIC
  aics <- rep(Inf, length(lambda))
  for (k in 1:length(lambda)) {
    # Sigma using k principle components on obsGrid
    Sigma_Y <- phiObs[,1:k,drop=FALSE] %*% diag(lambda[1:k], nrow = k) %*% t(phiObs[,1:k,drop = FALSE]) + diag(sigma2, nrow(phiObs))
    # get subset and index for each i
    ret <- lapply(t, function(tvec) {
      ind <- match(tvec, obsGrid)
      if (sum(is.na(ind)) != 0) {
        stop('Time point not found in obsGrid. ')
      }
      return(list(muVec=muObs[ind], Sigma_Yi=Sigma_Y[ind,ind,drop=FALSE]))
    })
    if (method == "Li") {
      sigma2k <- mapply(function(yVec, reti) {
        sum((sigma2*solve(reti$Sigma_Yi) %*% (yVec-reti$muVec))^2)
      }, y, ret)
      aics[k] <- N*log(sum(sigma2k)/N) + N + 2*n*k 
    } else {
      sigma2k <- mapply(function(yVec, reti) {
        solve(reti$Sigma_Yi, yVec-reti$muVec) %*% (yVec-reti$muVec) + log(det(reti$Sigma_Yi))
      }, y, ret)
      aics[k] <- sum(sigma2k) + 2*k
    }
    if (k>1 && aics[k]>aics[k-1]) { # stop decreasing 
      K <- k-1;
      criterion <- aics[k-1];
      break;
    } else if (k == length(lambda)) { # select the smallest AIC
      K <- which.min(aics)
      criterion <- min(aics)
    }
  }

  ## fitted covariance matrix on obsGrid and workGrid
  fittedCovObs <- phiObs[,1:K,drop=FALSE] %*% diag(lambda[1:K],K) %*% t(phiObs[,1:K,drop=FALSE])
  fittedCov <- phi[,1:K,drop=FALSE] %*% diag(lambda[1:K],K) %*% t(phi[,1:K,drop=FALSE])
  
  return(list(lambda = lambda, phi = phi, phiObs = phiObs[,1:K,drop=FALSE], fittedCovObs = fittedCovObs,
              fittedCov = fittedCov, kChoosen = K, criterion = criterion, cumFVE = FVE))
}


# test Yao's AIC: bad when measurement error is small
# Ly <- split(W_s, as.factor(rep(1:n, each = m)))
# Lt <- lapply(1:n, function(x) {s[x,]})
# pca <- FPCA(Ly, Lt, optns = list(dataType = "Sparse", FVEthreshold = 1, nRegGrid=grid.num))
# muObs <- FPCAIC(Ly, Lt, optns = list(dataType = "Sparse", FVEthreshold = 1, nRegGrid=grid.num),
#                 spec = list("AIC", FALSE))$muObs
# phiObs <- ConvertSupport(pca$workGrid, pca$obsGrid, phi=pca$phi)
# obsGrid <- sort(unique( c(unlist(Lt))))
# aics <- rep(0, length(pca$selectK))
# for (k in 1:length(pca$lambda)) {
#   phixi <- phiObs[,1:k,drop=FALSE] %*% t(pca$xiEst[,1:k,drop=FALSE])
#   # get subset and index for each i
#   ret <- lapply(Lt, function(tvec) {
#     ind <- match(tvec, obsGrid)
#     if (sum(is.na(ind)) != 0) {
#       stop('Time point not found in obsGrid. ')
#     }
#     return(list(muVec=muObs[ind], phixi = phixi[ind]))
#   })
#   sigma2k <- mapply(function(yVec, reti) {
#     sum((yVec-reti$muVec-reti$phixi)^2)
#   }, Ly, ret)
#   aics[k] <- sum(sigma2k)/pca$sigma2 + 2*k
# }

