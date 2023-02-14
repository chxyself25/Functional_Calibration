# simulation for time-varying case, setting I and setting II
## compare with simple last observation carried forward and Cao's paper
#source("../Funcs/sim_funcs.R")
#source("../Funcs/Asyncreg_X1.R")
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
environment(GetCrCovYX) <- asNamespace("fdapace")

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
bb = NULL

## simulation 200 times
# res <- foreach(sim = 1:200, .combine = "rbind") %dopar% {
res <- NULL
for (batch in 1:20) {
cat(batch, "\n")
res_list <- furrr::future_map(1:10, ~{
  ## univariate representation simulation
  data.list <- sim_dataxt1(n=n, m=m, sd.pcs = sd.pcs, sd.e = 0)
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
  Lw <- split(W_s, as.factor(rep(1:n, each = m)))
  Ls <- lapply(1:n, function(x) {s[[1]][x,]})
  #Ly <- split(Y, as.factor(rep(1:n, each = m)))
  t0 <- Sys.time()
  fm.res <- tvFMlassoX1(Y, t, Lw, Ls, optns, convert = TRUE, bw = bb)
  fm.time <- as.numeric(Sys.time()-t0, units = "mins")
  # simple last observation carried forward, same result with the package
  # locf.res <- tvLOCFlassoX1(Y, t, Lw, Ls, grid.len = optns$nRegGrid, bw = bb)
  # kernel last observation carried forward
  t0 <- Sys.time()
  czf.res <- tvCZFlassoX1(Y, t, W_s, s[[1]], grid.len = optns$nRegGrid, bw = bb) 
  czf.time <- as.numeric(Sys.time()-t0, units = "mins")
  #data.x <- data.frame(ID = rep(1:n, each = m), s = c(t(s[[1]])), W = W_s)
  #data.y <- data.frame(ID = rep(1:n, each = m), t = c(t(t)), Y = Y)
  #ti.res <- asynchTD(data.x, data.y, times = fm.res$time, kType = "epan", lType = "identity", bw = 1)
  #b.ti <- ti.res$betaHat
  # time-varying coefficient function
  t0 <- Sys.time()
  f.res <- tvFunctionlassoX1(Y, t, Lw, Ls, optns, optns, cr.bw = c(bb, bb))
  f.time <- as.numeric(Sys.time()-t0, units = "mins")
  # oracle estimate
  t0 <- Sys.time()
  o.res <- tvOracle(Y, t, X_t, grid.len = optns$nRegGrid, bw = bb)
  o.time <- as.numeric(Sys.time()-t0, units = "mins")
  # summarize results
  resdf <- NULL
  #for (i in c("mae", "mse", "maeint", "mseint")) {
  #  dfm <- data.frame(value = c(fm.res[[i]], locf.res[[i]], czf.res[[i]], f.res[[i]]), metric = rep(paste(i, c(0,1), sep = ""),4), 
  #                    method = rep(c("fImpute", "LOCF", "CZF", "FUNC"), each = 2), id = sim)
  #  resdf <- rbind(resdf, dfm)
  #}
  for (i in c("fm", "czf", "f", "o")) {
    resi <- get(paste0(i, ".res"))
    timi <- get(paste0(i, ".time"))
    resdf <- rbind(resdf, data.frame(id = .x+10*(batch-1), method = i, 
                                     beta0 = resi$coef[,1], beta1 = resi$coef[,2], time = resi$time, 
                                     bw = resi$bw, run_time = timi))
  }
  resdf
}, .options = furrr_options(seed = TRUE))

res <- rbind(res, do.call("rbind", res_list))
saveRDS(res, file = "./onext1_sim_nind_noerror_200.rds")
}

print(warnings())
