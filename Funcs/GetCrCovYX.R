

GetCrCovYX <- function(bw1 = NULL, bw2 = NULL, Ly1, Lt1 = NULL, Ymu1 = NULL, Ly2, Lt2 = NULL, Ymu2 = NULL, 
                       useGAM = FALSE, rmDiag=FALSE, kern='gauss', bwRoutine = 'l-bfgs-b', grid1=NULL, grid2=NULL) {
  
  if (kern != 'gauss' && (is.null(bw1) || is.null(bw2))) {
    stop('Cannot select bandwidth for non-Gaussian kernel')
  }
  
  if( !(bwRoutine %in% c('grid-search', 'bobyqa', 'l-bfgs-b') ) ){ 
    stop("'bwRoutine' argument is invalid.")
  }
  
  # If only Ly1 and Ly2 are available assume DENSE data
  if( is.matrix(Ly1) && is.null(Lt1) && is.null(Ymu1) && is.matrix(Ly2) && is.null(Lt2) && is.null(Ymu2)){
    rawCC <- GetRawCrCovFuncFunc(Ly1 = Ly1, Ly2 = Ly2)
    return ( list(smoothedCC = NULL, rawCC = rawCC, bw = NULL, score = NULL) )  
  }
  
  # Otherwise assume you have SPARSE data
  if( is.null(Ymu1) ||   is.null(Ymu2)){
    stop("Both functional means must be provided.")   
  }  
  
  # Get the Raw Cross-covariance    
  rawCC = GetRawCrCovFuncFunc(Ly1 = Ly1, Lt1 = Lt1, Ymu1 = Ymu1, Ly2 = Ly2, Lt2 = Lt2, Ymu2 = Ymu2)
  # Use a heuristic to decide when to bin
  if (!useGAM && sum(duplicated(rawCC$tpairn)) >= 0.2 * length(rawCC$rawCCov)) {
    # message('Binning rawCC')
    rawCC <- BinRawCov(rawCC)
  } else {
    # rename the fields to be the same as the binned rawCC
    names(rawCC)[names(rawCC) == 'rawCCov'] <- 'meanVals'
    names(rawCC)[names(rawCC) == 'tpairn'] <- 'tPairs'
    rawCC$count <- rep(1, length(rawCC$meanVals))
  }
  
  if (rmDiag) {
    diagInd <- rawCC$tPairs[, 1] == rawCC$tPairs[, 2]
    rawCC$tDiag <- rawCC$tPairs[diagInd, , drop=FALSE]
    rawCC$diagMeans <- rawCC$meanVals[diagInd]
    rawCC$diagCount <- rawCC$count[diagInd]
    rawCC$diagRSS <- rawCC$RSS[diagInd]
    rawCC$tPairs <- rawCC$tPairs[!diagInd, , drop=FALSE]
    rawCC$meanVals <- rawCC$meanVals[!diagInd]
    rawCC$count <- rawCC$count[!diagInd]
    rawCC$RSS <- rawCC$RSS[!diagInd]
  }
  
  # Calculate the observation and the working grids
  ulLt1 = unlist(Lt1);             ulLt2 = unlist(Lt2)
  obsGrid1 = sort(unique(ulLt1));  obsGrid2 = sort(unique(ulLt2))
  
  if (is.null(grid1) || is.null(grid2)) {
    workGrid1 = seq(min(obsGrid1), max(obsGrid1), length.out = 60)
    workGrid2 = seq(min(obsGrid2), max(obsGrid2), length.out = 60)
    grid.len <- 60
  } else {
    workGrid1 <- grid1
    workGrid2 <- grid2
    grid.len <- length(grid1)
  }
  workGrid12 = matrix(c(workGrid1, workGrid2),ncol= 2)
  
  if (useGAM == TRUE){ 
    # stop('Cannot be used on binned rawCC') # we cannot use binning for useGAM
    Qdata = data.frame(x =  rawCC$tPairs[,1], y = rawCC$tPairs[,2], z = rawCC$meanVals, group = rawCC$IDs  )
    # I comparsed with 're', ds' and 'gp' too, and 'tp' seems to have a slight edge for what we want
    # myGAM = mgcv::gamm( z ~ s(x,y, bs =c('tp','tp')), random=list(group=~1) , data= Qdata)$gam
    myGAM = mgcv::gam( z ~ s(x,y, bs =c('tp','tp')), data= Qdata)
    estPoints = data.frame( x= rep(workGrid1, times=grid.len), y= rep(workGrid2, each =grid.len), group = rep(3,grid.len*grid.len )) 
    smoothedCC = matrix(predict(myGAM, estPoints), grid.len)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = NULL, score = NULL) )  
  }
  
  # If the bandwidth is known already smooth the raw CrCov
  if( is.numeric(bw1) &&  is.numeric(bw2)){
    smoothedCC <- smoothRCC2D(rcov =rawCC, bw1, bw2, workGrid1, workGrid2,
                              kern=kern)    
    # potentially incorrect GCV score if the kernel is non-Gaussian
    score = GCVgauss2D(smoothedCC = smoothedCC, smoothGrid = workGrid12, 
                       rawCC = rawCC$meanVals, rawGrid = rawCC$tPairs, 
                       bw1 = bw1, bw2 = bw2)                      
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw =  c(bw1, bw2), score = score, smoothGrid = workGrid12 ) )
    
    # If the bandwidths are unknown use GCV to take find it
  } else {
    # Construct candidate bw's
    bwCandidates <- getBWidths(Lt1, Lt2)
    minGcvScores = Inf
    
    if(bwRoutine == 'grid-search'){
      # Find their associated GCV scores 
      gcvScores = rep(Inf, nrow(bwCandidates)) 
      for (i in 1:length(bwCandidates)){
        gcvScores[i] = theCostFunc(bwCandidates[i,], rawCC, workGrid1, workGrid2, kern, workGrid12)  
      } 
      # Pick the one with the smallest score 
      minimumScore =  min(gcvScores, na.rm=TRUE)
      if( minimumScore == Inf) {
        stop("It seems that the minimum GCV score equals infinity. Stopping 'GetCrCovYX'")
      }
      bInd = which(gcvScores == minimumScore);
      bOpt1 = max(bwCandidates[bInd,1]);
      bOpt2 = max(bwCandidates[bInd,2]); 
      minGcvScores = minimumScore
    } else {
      
      bwRanges = apply(bwCandidates,2, range)
      upperB = bwRanges[2,]
      lowerB = bwRanges[1,]
      
      if( !is.element('minqa', installed.packages()[,1]) && bwRoutine == 'bobyqa'){
        warning("Cannot use 'minqa::bobyqa' to find the optimal bandwidths. 'minqa' is not installed. We will do an 'L-BFGS-B' search.")
        bwRoutine == 'l-bfgs-b'
      }
      
      if( bwRoutine == 'l-bfgs-b' ){
        theSols = optim(fn = theCostFunc, par = upperB*0.95, # Starting value that "is safe"
                        upper = upperB, lower = lowerB, method ='L-BFGS-B', control = list(maxit = grid.len),
                        rawCC = rawCC, workGrid1 = workGrid1, workGrid2 = workGrid2, kern = kern, workGrid12 = workGrid12) 
        minGcvScores = theSols$value
      } else { # when BOBYQA is available
        theSols = minqa::bobyqa(fn = theCostFunc, par = upperB*0.95, # Starting value that "is safe"
                                upper = upperB, lower = lowerB, control = list(maxfun = 41),
                                rawCC, workGrid1, workGrid2, kern, workGrid12) 
        minGcvScores = theSols$fval
      }
      
      bOpt1 = theSols$par[1]
      bOpt2 = theSols$par[2]
      
      if( bOpt1 > 0.75 * upperB[1] && bOpt2 > 0.75 * upperB[2] ){
        warning('It seems that the bandwidth selected by the solvers is somewhat large. Maybe you are in local minima.')
      }
    }
    smoothedCC <- smoothRCC2D(rcov=rawCC, bw1 =bOpt1, bw2 =bOpt2, workGrid1, workGrid2, kern=kern)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = c(bOpt1, bOpt2), smoothGrid = workGrid12, 
                  score = minGcvScores) )
  }  
}
