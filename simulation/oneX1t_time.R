#### timing for a typical run of time-varying coefficient model #########
library(tvReg) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(MASS)
library(AsynchLong) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(caTools) #, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
#library(doParallel)
#registerDoParallel(cores = 10)
library(future)
library(furrr)
future::plan(multicore, workers = 10)
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../Funcs/")
environment(kernelAuto_wt) <- asNamespace("AsynchLong")
environment(asynchTD) <- asNamespace("AsynchLong")

## specification
# mean function
mean_fun <- function(tt) {
  tt + sin(tt)
}
# define eigen functions: k = 1,2,3
psi <- function(t, k) {
  (1/sqrt(5))*sin(pi*k*t/10)
}

# error variance structure
evar_fun <- function(tt, ind = FALSE) {
  if (ind) {
    covs <- diag(1.5, nrow = length(tt))
  } else {
    covs <- outer(tt, tt, function(t1,t2) {
      2^(-abs(t1-t2)/5)
    })}
  return(covs)
}
# time-dependent beta: including setting I and II
beta_fun <- function(tt, intercept = FALSE) {
  if (intercept) {
    betat <- 0.2*tt + 0.5 # setting I
    #betat <- tt^(0.5) # setting II
  } else {
    betat <- sin(pi*tt/10) # setting I 
    #betat <- sin(pi*tt/5) # setting II
  }
  return(betat)
}

n <- 200; m <- 5
sd.pcs <- c(2, sqrt(2), 1)
optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "GCV", methodBwCov = "GCV")
#bb <- 1.034357
#bb <- 2.355484

set.seed(1000)
data.list <- sim_dataxt1(n=n, m=m, sd.pcs = sd.pcs, sd.e = 1)
#saveRDS(data.list, file = paste0("./sim_data/onext2_sim_", sim, ".rds"))
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
bb = NULL
Lw <- split(W_s, as.factor(rep(1:n, each = m)))
Ls <- lapply(1:n, function(x) {s[[1]][x,]})
#Ly <- split(Y, as.factor(rep(1:n, each = m)))
t0 <- Sys.time()
fm.res <- tvFMlassoX1(Y, t, Lw, Ls, optns, convert = TRUE, bw = bb)
fm.time <- as.numeric(Sys.time()-t0, units = "mins")
mubw = fm.res$pca.mubw; covbw = fm.res$pca.covbw; bb = fm.res$bw
### bootstrap for FCAR
t0 <- Sys.time()
b.bts_list <- furrr::future_map(1:10, ~{
  bts_res = NULL
  for (bts_ in 1:50) {
    s.idx <- sample(1:n, size = n, replace = TRUE)
    Lwb <- Lw[s.idx]
    Lsb <- Ls[s.idx]
    tb <- t[s.idx,]
    Yb <- c(matrix(Y, ncol = n, byrow = FALSE)[,s.idx])
    fm.b <- tvFMlassoX1(Yb, tb, Lwb, Lsb, optns = list(dataType = "Sparse", nRegGrid = 60, userBwMu = mubw, userBwCov = covbw), 
                        convert = TRUE, bw = bb)
    #c(fm.b$beta) 
    bts_res <- rbind(bts_res, data.frame(beta0 = fm.b$coef[,1], beta1 = fm.b$coef[,2], time = fm.b$time, 
                                         iter = bts_ + 50*(.x-1)))
  }
  bts_res
}, .options = furrr_options(seed = TRUE))
b.bts <- do.call("rbind", b.bts_list); rm(b.bts_list)
sd_time.fm <- as.numeric(Sys.time() - t0, units = "mins")

# kernel last observation carried forward
bb = NULL
#czf.res <- tvCZFlassoX1(Y, t, W_s, s[[1]], grid.len = optns$nRegGrid, bw = bb) 
n <- nrow(t)
m <- ncol(t)
## rescale s and t
ss <- s[[1]]/10; tt <- t/10
data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(ss)), W = W_s)
data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(tt)), Y = Y)
#timepts <- sort(c(t))
timepts <- seq(max(min(tt), min(ss)), min(max(tt), max(ss)), length.out = optns$nRegGrid)
#ti.res <- foreach(x = timepts, .combine = "rbind") %do% {
t0 <- Sys.time()
ti.res <- lapply(timepts, function(x) {
  timodel <- tryCatch(asynchTD(data.x, data.y, times = x, kType = "epan", lType = "identity",
                               bw = bb, verbose = FALSE), error = function(e) {NA})
  if (any(is.na(timodel))) {
    rep(NA, 4)
  } else {
    c(c(timodel$betaHat), timodel$beta_time, timodel$sd_time) 
  }
})
ti.res <- do.call("rbind", ti.res)
czf.time <- as.numeric(Sys.time()-t0, units = "mins")
beta_time <- ti.res[,3]
beta_time[is.na(beta_time)] <- mean(ti.res[,3], na.rm = TRUE)
sum(beta_time)
sd_time <- ti.res[,4]
sd_time[is.na(sd_time)] <- mean(ti.res[,4], na.rm = TRUE)
sum(sd_time)
