#### get optimal FPCA bandwidth #####
library(MASS)
library(Matrix)
library(fdapace)
library(dplyr)
library(doParallel)
registerDoParallel(cores = 20)
sourceDir <- function(path, trace = FALSE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../Funcs/")

############ load the data  #############
cat("loading Xvars data file", "\n")
Xvars <- readRDS("./Xvars_selected.rds")

######### regression analysis using multivariate functional calibration ###############
optns <- list(dataType = "Sparse", nRegGrid = 100, methodBwMu = "GCV", methodBwCov = "GCV")
types <- c("CVR", "PHY")
vars <- c("TRIGRES", "BMI")

arg1 <- paste0("!is.na(", vars[1], ") & ", vars[1], ">0")
arg2 <- paste0("!is.na(", vars[2], ") & ", vars[2], ">0")
dat1 <- Xvars %>% filter(eval(parse(text = paste0("!is.na(", types[1], "DAY) & ", types[1], "DAY >= 0")))) %>% 
  filter(eval(parse(text = arg1))) %>% as.data.frame
dat2 <- Xvars %>% filter(eval(parse(text = paste0("!is.na(", types[2], "DAY) & ", types[2], "DAY >= 0")))) %>%
  filter(eval(parse(text = arg2))) %>% as.data.frame
dat1.list <- split(dat1, f = as.factor(dat1$ID)); nms1 <- names(dat1.list)
dat2.list <- split(dat2, f = as.factor(dat2$ID)); nms2 <- names(dat2.list)

dx <- length(vars)
W_s.list <- list()
Lsi <- lapply(dat1.list, function(x) {x[[paste0(types[1], "DAY")]]/365})
Lwi <- lapply(dat1.list, '[[', vars[1])
W_s.list[[1]] <- remove_duplicates(Lsi, Lwi)

Lsi <- lapply(dat2.list, function(x) {x[[paste0(types[2], "DAY")]]/365})
Lwi <- lapply(dat2.list, '[[', vars[2])
W_s.list[[2]] <- remove_duplicates(Lsi, Lwi)
## univariate FPCA and covariance calculation
cat("univariate FPCA:", "\n")
pca.list = list()
for (v in 1:length(W_s.list)) {
  cat(vars[v], "\n")
  Ly <- W_s.list[[v]]$Ly
  Lt <- W_s.list[[v]]$Lt
  pca <- FPCAIC(Ly, Lt, optns, spec = list("Li", FALSE))
  pca.list[[v]] <- pca
}
bw.res <- lapply(pca.list, function(x) {x[c("bwMu", "bwCov")]})
saveRDS(bw.res, "./bw_results_GCV.rds")

