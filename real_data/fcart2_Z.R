#### regression analysis of 2 covariates, adjusting age and income ##########
#### using multivariate functional calibration ##############################
library(MASS)
library(Matrix)
library(fdapace)
library(tvReg)
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
environment(GetCrCovYX) <- asNamespace("fdapace")

############ load the data  #############
cat("loading Xvars data file", "\n")
Xvars <- readRDS("./Xvars_selected.rds")
IDs <- unique(Xvars$ID)
#dat0 <- Xvars %>% filter(!is.na(COGDAY) & COGDAY >= 0) %>% 
#  filter(!is.na(GLBSCORE)) %>% as.data.frame
#Lt <- lapply(split(dat0, f = as.factor(dat0$ID)), '[[', 'COGDAY')
#Ly <- lapply(split(dat0, f = as.factor(dat0$ID)), '[[', 'GLBSCORE')
dat0 <- Xvars %>% filter(!is.na(HRMDAY) & HRMDAY >= 0) %>% 
  filter(!is.na(FSH)) %>% as.data.frame
Lt <- lapply(split(dat0, f = as.factor(dat0$ID)), function(x) {x[['HRMDAY']]/365})
Ly <- lapply(split(dat0, f = as.factor(dat0$ID)), '[[', 'FSH')
tylist <- remove_duplicates(Lt, Ly)
Lt <- tylist$Lt; Ly <- tylist$Ly
#nmt <- names(Lt)
Zvars <- readRDS("./Zvars_selected.rds") %>% filter(!is.na(income_level) & !is.na(AGE)) %>% 
  mutate(SWANID = as.character(SWANID)) %>% as.data.frame
nmt <- intersect(names(Ly), Zvars$SWANID)
Lt <- Lt[nmt]; Ly <- Ly[nmt]
Zvars$income_level <- factor(Zvars$income_level, levels = c("low", "high"))
#Zvars <- Zvars %>% filter(SWANID %in% nmt) %>% 
#  mutate(race = if_else(RACE %in% c("2", "3"), "2_3", RACE)) %>% as.data.frame
#Zvars$race <- as.factor(Zvars$race)
Zvars_exp <- Zvars[match(rep(nmt, times = sapply(Ly, length)), Zvars$SWANID),]
#z_mat <- model.matrix(~0+race, Zvars_exp)
z_mat <- model.matrix(~1+AGE+income_level, Zvars_exp)
colnames(z_mat) <- c("intercept", "age", "income_high")

######### regression analysis using multivariate functional calibration ###############
pca.list <- readRDS("./bw_results_GCV.rds")
#optns <- list(dataType = "Sparse", nRegGrid = 100, userBwMu = 350/365, userBwCov = 500/365)
#bbs <- c(600, 600)/365
#optns <- list(dataType = "Sparse", nRegGrid = 60, methodBwMu = "GCV", methodBwCov = "GCV")
#types <- c("HRM", "CVR")
types <- c("CVR", "PHY")
res_coef <- list()
res_pca <- list()
res_bs <- list()
for (h in c("TRIGRES")) {
  cat("doing ", h, "\n")
  vars <- c(h, "BMI")
  # make sure observation day is positive
  arg1 <- paste0("!is.na(", vars[1], ") & ", vars[1], ">0")
  arg2 <- paste0("!is.na(", vars[2], ") & ", vars[2], ">0")
  dat1 <- Xvars %>% filter(eval(parse(text = paste0("!is.na(", types[1], "DAY) & ", types[1], "DAY >= 0")))) %>% 
    filter(eval(parse(text = arg1))) %>% as.data.frame
  dat2 <- Xvars %>% filter(eval(parse(text = paste0("!is.na(", types[2], "DAY) & ", types[2], "DAY >= 0")))) %>%
    filter(eval(parse(text = arg2))) %>% as.data.frame
  #if (h == "DHAS") {
  #  res_coef <- readRDS("./fcart_mr_GLUCOSE_HRM_fix_sd.rds")
  #  res_pca <- readRDS("./fcart_mr_GLUCOSE_HRM_pca_fix_sd.rds")
  #  coefh <- res_coef[[h]]
  #} else {
  dat1.list <- split(dat1, f = as.factor(dat1$ID)); nms1 <- names(dat1.list); subject_n1 = length(dat1.list)
  dat2.list <- split(dat2, f = as.factor(dat2$ID)); nms2 <- names(dat2.list); subject_n2 = length(dat2.list)
  optns1 <- list(dataType = "Sparse", nRegGrid = 100, kernel = "gauss",
                 userBwMu = unname(pca.list[[1]]$bwMu)*subject_n1^(-1/10), userBwCov = pca.list[[1]]$bwCov*subject_n1^(-1/10))
  optns2 <- list(dataType = "Sparse", nRegGrid = 100, kernel = "gauss",
                 userBwMu = unname(pca.list[[2]]$bwMu)*subject_n2^(-1/10), userBwCov = pca.list[[2]]$bwCov*subject_n2^(-1/10))
  optns <- list(optns1, optns2)
  bbs <- rep(max(optns1$userBwCov,optns2$userBwCov)*(length(Ly))^(1/40), 2)
  resh <- fcart_multi(dat1.list, dat2.list, z_mat, Lt, Ly, types, vars, optns, bws = bbs)
  res_coef[[h]] <- coefh <- resh[[1]]
  res_pca[[h]] <- resh[[2]]
  saveRDS(res_coef, file = "./fcart_mr_FSH_BMI_INCOME_AGE_fix_us_Li.rds")
  saveRDS(res_pca, file = "./fcart_mr_FSH_BMI_INCOME_AGE_pca_fix_us_Li.rds") 
  #}
  ## do bootstrap 
  #optns$userBwMu <- resh[[2]]$pca.list$bwMu
  #optns$userBwCov <- resh[[2]]$pca.list$bwCov
  if (TRUE) {
    bws = c(coefh$bw[coefh$method == "FCART"][1], coefh$bw[coefh$method == "LOCF"][1])
    cat("start bootstrap:", "\n")
    resh <- NULL
    for (b in 1:200) {
      cat("bootstrap at " ,b, "\n") 
      #resh <- foreach(b = 1:200, .combine = "rbind") %dopar% {
      #s.idx <- sample(1:n, n, replace = TRUE)
      s.id <- sample(IDs, length(IDs), replace = TRUE)
      Lyb <- Ly[s.id[s.id %in% nmt]]
      Ltb <- Lt[s.id[s.id %in% nmt]]
      Zvars_expb <- Zvars[match(rep(names(Lyb), times = sapply(Lyb, length)), Zvars$SWANID),]
      z_matb <- model.matrix(~1+AGE+income_level, Zvars_expb)
      colnames(z_matb) <- c("intercept", "age", "income_high")
      datb1.list <- dat1.list[s.id[s.id %in% nms1]]
      datb2.list <- dat2.list[s.id[s.id %in% nms2]]
      #saveRDS(list(datb1.list, datb2.list, z_matb, Ltb, Lyb), file = "./temp_bootstrap.rds")
      resb <- fcart_multi(datb1.list, datb2.list, z_matb, Ltb, Lyb, types, vars, optns, bws = bws)
      resb[[1]]$bts <- b
      resh <- rbind(resh, resb[[1]])
    }
    res_bs[[h]] <- resh
    saveRDS(res_bs, file = "./fcart_mr_FSH_BMI_INCOME_AGE_bs200_fix_us_Li.rds")
  }
}
saveRDS(res_coef, file = "./fcart_mr_FSH_BMI_INCOME_AGE_fix_us_Li.rds")
saveRDS(res_pca, file = "./fcart_mr_FSH_BMI_INCOME_AGE_pca_fix_us_Li.rds")
saveRDS(res_bs, file = "./fcart_mr_FSH_BMI_INCOME_AGE_bs200_fix_us_Li.rds")


