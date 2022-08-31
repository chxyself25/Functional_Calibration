# mean function part of FPCA, used for centering univariate functional data

FPCAMu <- function(Ly, Lt, optns = list()) {

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

  ret <- list(obsGrid = obsGrid, muObs = muObs, workGrid = workGrid, muWork = muWork, resid = resid, bwMu = smcObj$bw_mu)
  return(ret)
}
