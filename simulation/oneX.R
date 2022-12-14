## simulation for time-invariant case setting I
## compare with simple last observation carried forward and Cao's paper
#source("../Funcs/sim_funcs.R")
#source("../Funcs/Asyncreg_X1.R")
library(fdapace) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(MASS)
library(AsynchLong) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(caTools) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
#library(doParallel)
#registerDoParallel(cores = 10)
library(future)
library(furrr)
future::plan(multisession, workers = 10)
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../Funcs/")
environment(kernelAuto_wt) <- asNamespace("AsynchLong")
environment(asynchTI) <- asNamespace("AsynchLong")


## specification
# mean function
mean_fun <- function(tt) {
  tt + sin(tt)
}
# define eigen functions: k = 1,2,3
psi <- function(t, k) {
  (1/sqrt(5))*sin(pi*k*t/10)
}

# error variance structure: ind = TRUE for independent error, ind = FALSE for non-independent error
evar_fun <- function(tt, ind = FALSE) {
  if (ind) {
    covs <- diag(1.5, nrow = length(tt))
  } else {
  covs <- outer(tt, tt, function(t1,t2) {
    2^(-abs(t1-t2)/5)
  })}
  return(covs)
}

b <- c(1, 2); n <- 200; m <- 5
sd.pcs <- c(2, sqrt(2), 1)
optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "GCV", methodBwCov = "GCV")

## test for simulation function
##data.list <- sim_datax(n=n, m=m, sd.pcs = sd.pcs, sd.e = 1, dx = 1, beta = b, method = "fix")

## simulation 200 times
res <- NULL
for (i in 1:200) {
  cat("Doing the simulation ", i, "th time", "\n")
  ## univariate representation simulation
  data.list <- sim_datax(n=n, m=m, sd.pcs = sd.pcs, sd.e = 0, dx = 1, beta = b, method = "fix")
  # observed X (on S domain)
  W_s <- data.list$Ws
  # true X on t domain
  X_t <- data.list$Xt
  # two observed sampling time points
  s <- data.list$s
  t <- data.list$t
  # observed y response
  Y <- data.list$Y
  
  # functional match 
  Lw <- split(W_s, as.factor(rep(1:n, each = m)))
  Ls <- lapply(1:n, function(x) {s[[1]][x,]})
  t0 <- Sys.time()
  fm.res <- FMlassoX1(Y, t, Lw, Ls, optns = optns)
  est_time.fm <- as.numeric(Sys.time()-t0, units = "mins")
  b.fm <- fm.res$beta
  mubw <- fm.res$pca.mubw; covbw <- fm.res$pca.covbw
  # bootstrap sampling for beta standard error estimation
  t0 = Sys.time()
  #b.bts <- foreach(bts = 1:500, .combine = "rbind") %dopar% {
  b.bts_list <- furrr::future_map(1:10, ~{
    bts_res = c()
    for (bts_ in 1:50) {
      s.idx <- sample(1:n, size = n, replace = TRUE)
      Lwb <- Lw[s.idx]
      Lsb <- Ls[s.idx]
      tb <- t[s.idx,]
      Yb <- c(matrix(Y, ncol = n, byrow = FALSE)[,s.idx])
      fm.b <- FMlassoX1(Yb, tb, Lwb, Lsb, optns = list(dataType = "Sparse", nRegGrid = 60, userBwMu = mubw, userBwCov = covbw))
      #c(fm.b$beta) 
      bts_res <- rbind(bts_res, c(fm.b$beta))
    }
    bts_res
  }, .options = furrr_options(seed = TRUE))
  b.bts <- do.call("rbind", b.bts_list); rm(b.bts_list)
  sd_time.fm <- as.numeric(Sys.time() - t0, units = "mins")
  # simple last observation carried forward, same result with the package
  # locf.res <- LOCFlassoX(Y, t, W_s, s)
  # b.locf <- locf.res$beta
  # kernel last observation carried forward
  data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(s[[1]])), W = W_s)
  data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(t)), Y = Y)
  ti.res <- asynchTI(data.x, data.y)
  b.ti <- ti.res$betaHat
  # summarize results
  # betas <- cbind(b.fm, b.locf, as.matrix(b.ti))
  # resi <- data.frame(beta0 = betas[1,], beta1 = betas[2,], method = c("fm", "locf", "ti"),
  #            pca.K = c(fm.res$pca.K, NA, NA), pca.sigma2 = c(fm.res$pca.sigma2, NA, NA), 
  #            pca.rho = c(fm.res$pca.rho, NA, NA), naive.sigma20 = c(fm.res$beta.sigma2[1], locf.res$beta.sigma2[1], NA),
  #            naive.sigma21 = c(fm.res$beta.sigma2[2], locf.res$beta.sigma2[2], NA), std0 = c(sd(b.bts[,1]), NA, ti.res$stdErr[1]),
  #            std1 = c(sd(b.bts[,2]), NA, ti.res$stdErr[2]))
  betas <- cbind(b.fm, as.matrix(b.ti))
  resi <- data.frame(beta0 = betas[1,], beta1 = betas[2,], method = c("fm", "kw"),
                     pca.K = c(fm.res$pca.K, NA), pca.sigma2 = c(fm.res$pca.sigma2, NA), 
                     pca.rho = c(fm.res$pca.rho, NA), naive.sigma20 = c(fm.res$beta.sigma2[1], NA),
                     naive.sigma21 = c(fm.res$beta.sigma2[2], NA), std0 = c(sd(b.bts[,1]), ti.res$stdErr[1]),
                     std1 = c(sd(b.bts[,2]), ti.res$stdErr[2]), est_time = c(est_time.fm, ti.res$beta_time),
                     sd_time = c(sd_time.fm, ti.res$sd_time))
  res <- rbind(res, resi)
  saveRDS(res, file = paste0("./onex1_sim_nind_noerror_", n, ".rds"))
  cat("Have done the simulation ", i, "\n")
}
print(warnings())
#saveRDS(res, file = paste0("./onex1_sim_fix_", n, "_nind.rds")


