## FPCA for sparse data, but use Li's AIC to select number of principle components
## and use eigenfunctions and eigenvalues to reconstruct fitted covariance matrix
## a little simplified to save time, optns is the same, add one more argument AIC

FPCAIC = function(Ly, Lt, optns = list(), spec = list("AIC", TRUE)) {
  # Check the data validity for further analysis
  CheckData(Ly,Lt)
  
  # Force the data to be list of numeric members and handle NA's
  #Ly <- lapply(Ly, as.numeric) 
  #Lt <- lapply(Lt, as.numeric)
  #Lt <- lapply(Lt, signif, 14)
  #inputData <- list(Ly=Ly, Lt=Lt);
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  optns = SetOptions(Ly, Lt, optns);
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  CheckOptions(Lt, optns,numOfCurves)
  # m is too small, no binnning is necessary
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new
  # data
  # regGrid: the grid of time points for which the smoothed covariance
  # surface assumes values
  # cutRegGrid: truncated grid specified by optns$outPercent for the cov
  # functions
  obsGrid = sort(unique( c(unlist(Lt))));
  workGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(workGrid)
  
  ymat <- fdapace:::List2Mat(Ly, Lt)
  
  ## Mean function
  # If the user provided a mean function use it
  userMu <- optns$userMu
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    smcObj <- fdapace:::GetUserMeanCurve(optns, obsGrid, workGrid, buff)
    smcObj$muDense = ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
  } else if (optns$methodMuCovEst == 'smooth') { # smooth mean
    smcObj = fdapace:::GetSmoothedMeanCurve(Ly, Lt, obsGrid, workGrid, optns)
  }
  # mu: the smoothed mean curve evaluated at times 'obsGrid'
  muObs <- smcObj$mu
  resid <- lapply(1:numOfCurves, function(i) {
    Ly[[i]] - muObs[match(Lt[[i]], obsGrid)]
  })
  # mu: the smoothed mean curve evaluated at times 'workGrid'
  muWork <- smcObj$muDense

  ## Covariance function and sigma2
  if (!is.null(optns$userCov) && optns$methodMuCovEst != 'smooth') { 
    scsObj <- fdapace:::GetUserCov(optns, obsGrid, cutRegGrid, buff, ymat)
  } else if (optns$methodMuCovEst == 'smooth') {
    # smooth cov and/or sigma2
    scsObj = fdapace:::GetSmoothedCovarSurface(Ly, Lt, muObs, obsGrid, workGrid, optns,
                                     optns$useBinnedCov) 
  } 
  sigma2 <- scsObj[['sigma2']]
  # workGrid: possibly truncated version of the regGrid, regGrid if outpercent is 0, 1
  if (!all(workGrid == scsObj$outGrid)) {
    stop("Something is wrong, workGrid is truncated!")
  }
  
  # Get the results for the eigen-analysis: modified by using Li's AIC to choose number of principle components
  eigObj <- GetEigenAnalysisAIC(Ly, Lt, smoothCov=scsObj$smoothCov, workGrid, obsGrid, sigma2, muWork, muObs, 
                                spec, optns)
  CovObs <- eigObj$fittedCovObs; K <- eigObj$kChoosen
  CovObs.list <- lapply(1:numOfCurves, function(i) {
    idx <- match(Lt[[i]], obsGrid)
    CovObs[idx, idx, drop=FALSE]
  })
  phiObs <- eigObj$phiObs;
  phiObs.list <- lapply(1:numOfCurves, function(i) {
    phiObs[match(Lt[[i]], obsGrid),, drop=FALSE]
  })
  
  # Get scores
  if (optns$rho != 'no') { 
    if( is.null(optns$userRho) ){
      if( length(Ly) > 2048 ){
        randIndx <- sample( length(Ly), 2048)
        rho <- GetRho(Ly[randIndx], Lt[randIndx], optns, muObs, obsGrid, CovObs, eigObj$lambda[1:K], phiObs, sigma2)
      } else {
        rho <- GetRho(Ly, Lt, optns, muObs, obsGrid, CovObs, eigObj$lambda[1:K], phiObs, sigma2)
      }
    } else {
      rho = optns$userRho;
    }
    sigma2 <- rho
  }
  scoresObj <- GetCEScores(Ly, Lt, optns, muObs, obsGrid, CovObs, eigObj$lambda[1:K], phiObs, sigma2)
  
  # Make the return object: results are truncated at choosen K
  xiEst <- t(do.call(cbind, scoresObj[1, ])) 
  xiVar <- scoresObj[2, ] # variance of CE xi
  ret <- list(sigma2 = scsObj$sigma2, 
              lambda = eigObj$lambda, 
              phi = eigObj$phi,
              phiObs = phiObs.list,
              xiEst = xiEst, 
              xiVar = xiVar, 
              obsGrid = obsGrid, 
              workGrid = workGrid, 
              mu = muWork,
              #muObs = muObs, 
              smoothedCov = scsObj$smoothCov, 
              FVE = eigObj$cumFVE[eigObj$kChoosen], 
              cumFVE =  eigObj$cumFVE, 
              fittedCov = eigObj$fittedCov, 
              fittedCovObs = CovObs.list,
              optns = optns, 
              bwMu = smcObj$bw_mu, 
              bwCov = scsObj$bwCov,
              selectK = eigObj$kChoosen,
              criterion = eigObj$criterion,
              resid = resid)
  ret$rho <- if (optns$rho != 'no') rho else NULL
  ret$inputData <- inputData; # This will be potentially be NULL if `lean`
  class(ret) <- 'FPCA'
  
  return(ret); 
}
