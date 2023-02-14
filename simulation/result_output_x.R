########## summarize simulation results ################
## settings1: increasing trend mean, n = 200, independent/dependent error variance, no measurement error for covariates
## settings2, n = 200, independent/dependent error variance, no measurement error for covariates
library(dplyr)
# true beta and covariance of X in simulation
b0 <- 1; b1 <- 2; sd.pcs <- c(2, sqrt(2), 1)
select_types <- c("ind", "nind", "nind_noerror", "nind_m15")
fm.res <- NULL
for (set in 1:2) {
  for (type in select_types) {
    resi <- readRDS(paste0("./onex", set, "_sim_", type, "_200.rds"))
    fm.resi <- subset(resi, method == "fm") %>% 
      summarize(bias0 = mean(beta0-b0), est_std0 = sd(beta0-b0), 
                bias1 = mean(beta1-b1), est_std1 = sd(beta1-b1),
                naive0 = mean(sqrt(naive.sigma20)), naive.ci0 = mean(abs(beta0-b0) < 1.96*sqrt(naive.sigma20)), 
                naive1 = mean(sqrt(naive.sigma21)), naive.ci1 = mean(abs(beta1-b1) < 1.96*sqrt(naive.sigma21)),
                bts0 = mean(std0), bts.ci0 = mean(abs(beta0-b0) < 1.96*std0), 
                bts1 = mean(std1), bts.ci1 = mean(abs(beta1-b1) < 1.96*std1)) %>% as.data.frame
    fm.resi <- fm.resi %>% mutate_all(list(~round(., 3))) %>% mutate(setting = set, error_type = ifelse(type == "nind_m15_bwfix", "nind_m15", type)) %>% as.data.frame
    fm.res <- rbind(fm.res, fm.resi)
  }
}
out.names <- c("Bias", "SD", "Naive SE", "Naive CP", "Bootstrap SE", "Bootstrap CP")
select_metrics <- c("bias1", "est_std1", "naive1", "naive.ci1", "bts1", "bts.ci1")
#select_metrics <- c("bias0", "est_std0", "naive0", "naive.ci0", "bts0", "bts.ci0")
names(select_metrics) <- out.names
select_types <- c("ind", "nind", "nind_noerror", "nind_m15")
#select_types <- c("ind", "nind")
for (nm in out.names) {
  nm_value1 <- sapply(select_types, function(x) {
    fm.res[[select_metrics[nm]]][fm.res$setting == 1 & fm.res$error_type == x]
  })
  nm_value2 <- sapply(select_types, function(x) {
    fm.res[[select_metrics[nm]]][fm.res$setting == 2 & fm.res$error_type == x]
  })
  cat(paste(nm, paste(c(nm_value1, nm_value2), collapse = " & "), sep = " & "), "\\\\", "\n")
}


###### comparing with other methods: KW ######
select_types <- c("ind", "nind", "nind_noerror", "nind_m15")
#select_types <- c("ind", "nind")
out.names <- c("Bias", "SD", "SE", "CP")
select_metrics <- c("bias1", "est_std1", "mean.std1", "ci1")
#select_metrics <- c("bias0", "est_std0", "mean.std0", "ci0")
names(select_metrics) <- out.names
for (set in 1:2) {
  out_str <- ifelse(set == 1, "Setting I", "Setting \\II")
  set_tbl <- NULL
  for (type in select_types) {
    resi <- readRDS(paste0("./onex", set, "_sim_", type, "_200.rds"))
    resi$method <- gsub("ti", "kw", resi$method)
    #resi$method <- toupper(gsub("fm", "fcar", resi$method))
    all.resi <- subset(resi, method %in% c("fm", "kw")) %>% 
      group_by(method) %>% 
      summarize(bias0 = mean(beta0-b0), est_std0 = sd(beta0-b0),
                bias1 = mean(beta1-b1), est_std1 = sd(beta1-b1),
                mean.std0 = mean(std0), mean.std1 = mean(std1),
                ci0 = mean(abs(beta0-b0) < 1.96*std0), 
                ci1 = mean(abs(beta1-b1) < 1.96*std1)) %>% as.data.frame
    all.resi <- all.resi[match(c("fm", "kw"), all.resi$method), select_metrics]
    set_tbl <- cbind(set_tbl, round(t(all.resi), 3))
  }
  for (nm in out.names) {
    out_str <- paste0(paste(out_str, nm, sep = " & "), " & ", paste(set_tbl[select_metrics[nm],], collapse = " & "))
    cat(out_str, "\\\\", "\n")
    out_str <- ""
  }
}

#### computation time comparison ######
res_all <- NULL
for (set in 1:2) {
  resi <- readRDS(paste0("./onex", set, "_sim_nind_noerror_200.rds"))
  res_all <- rbind(res_all, resi %>% mutate(setting = set))
}
df <- res_all %>% 
  group_by(method) %>% 
  summarize(est_time_mean = mean(est_time),
            sd_time_mean = mean(sd_time)) %>% as.data.frame
df
