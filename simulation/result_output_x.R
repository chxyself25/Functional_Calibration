########## summarize simulation results ################
## settings1: increasing trend mean, n = 200, independent/dependent error variance
## settings2, n = 200, independent/dependent error variance
library(dplyr)
# true beta and covariance of X in simulation
b0 <- 1; b1 <- 2; sd.pcs <- c(2, sqrt(2), 1)
sigmaX <- sum(sd.pcs^2)
fm.res <- list()
for (set in 1:2) {
  for (ind in c("", "n")) {
    resi <- readRDS(paste0("./onex", set, "_sim_", ind, "ind_200.rds"))
    ## evaluation of functional match method
    fm.resi <- subset(resi, method == "fm") %>%
      summarize(SE = mean((beta0-b0)^2+(beta1-b1)^2), SE.std = sd((beta0-b0)^2+(beta1-b1)^2),
                PE = mean(sigmaX*(beta1-b1)^2), PE.std = sd(sigmaX*(beta1-b1)^2),
                bias0 = mean(beta0-b0), bias0.std = sd(beta0-b0), 
                bias1 = mean(beta1-b1), bias1.std = sd(beta1-b1),
                #beta.std0 = sd(beta0), beta.std1 = sd(beta1), 
                naive0 = mean(sqrt(naive.sigma20)), naive.ci0 = mean(abs(beta0-b0) < 1.96*sqrt(naive.sigma20)), 
                naive1 = mean(sqrt(naive.sigma21)), naive.ci1 = mean(abs(beta1-b1) < 1.96*sqrt(naive.sigma21)),
                bts0 = mean(std0), bts.ci0 = mean(abs(beta0-b0) < 1.96*std0), 
                bts1 = mean(std1), bts.ci1 = mean(abs(beta1-b1) < 1.96*std1)) %>% as.data.frame
    fm.resi <- apply(fm.resi, 2, round, 3)
    print(fm.resi)
    # output in format
    strsi <- c()
    for (j in seq(1, 8, by = 2)) {
      strsi <- c(strsi, paste0(fm.resi[j], "(", fm.resi[j+1], ")"))
    }
    names(strsi) <- names(fm.resi)[seq(1, 8, by = 2)]
    strsi <- c(strsi, fm.resi[9:16])
    fm.res[[paste0("set", set, ind)]] <- strsi
  }
}
tbl <- as.data.frame(fm.res)
# out.names <- c("SE", "PE", "$\\Delta \\beta_0$", "$\\Delta \\beta_x$", "$s.e.^*(\\wh \\beta_0)$", "$\\mathrm{CP}^*(\\wh \\beta_0)$",
#                "$s.e.^*(\\wh \\beta_x)$", "$\\mathrm{CP}^*(\\wh \\beta_x)$",
#                "$s.e.(\\wh \\beta_0)$", "$\\mathrm{CP}(\\wh \\beta_0)$", 
#                "$s.e.(\\wh \\beta_x)$", "$\\mathrm{CP}(\\wh \\beta_x)$")
# names(out.names) <- rownames(tbl)
# for (i in rownames(tbl)) {
#   cat(paste(out.names[i], paste(unlist(tbl[i,]), collapse = " & "), sep = " & "), "\\\\", "\n")
# }

###### comparing with other methods
out.names <- c("SE", "PE", "$\\Delta \\beta_0$", "$\\Delta \\beta_1$", 
               "$s.e.(\\wh \\beta_0)$", "$\\mathrm{CP}(\\wh \\beta_0)$",
               "$s.e.(\\wh \\beta_1)$", "$\\mathrm{CP}(\\wh \\beta_1)$")
for (set in 1:2) {
  cat("setting ", set, "\n")
  tbli <- list()
  for (ind in c("", "n")) {
    resi <- readRDS(paste0("./onex", set, "_sim_", ind, "ind_200.rds"))
    all.resi <- resi %>% 
      group_by(method) %>%
      summarize(SE = mean((beta0-b0)^2+(beta1-b1)^2), SE.std = sd((beta0-b0)^2+(beta1-b1)^2),
                PE = mean(sigmaX*(beta1-b1)^2), PE.std = sd(sigmaX*(beta1-b1)^2),
                bias0 = mean(beta0-b0), bias0.std = sd(beta0-b0),
                bias1 = mean(beta1-b1), bias1.std = sd(beta1-b1),
                mean.std0 = if_else(all(method == "locf"), mean(sqrt(naive.sigma20)), mean(std0)),
                ci0 = if_else(all(method == "locf"), mean(abs(beta0-b0) < 1.96*sqrt(naive.sigma20)), mean(abs(beta0-b0) < 1.96*std0)),
                mean.std1 = if_else(all(method == "locf"), mean(sqrt(naive.sigma21)), mean(std1)),
                ci1 = if_else(all(method == "locf"), mean(abs(beta1-b1) < 1.96*sqrt(naive.sigma21)), mean(abs(beta1-b1) < 1.96*std1))) %>% as.data.frame
    all.resi$method <- toupper(replace(as.character(all.resi$method), grep("ti", all.resi$method), "czf"))
    rownames(all.resi) <- all.resi$method; all.resi$method <- NULL
    all.resi <- round(all.resi, 3)
    for (m in c("FM", "CZF", "LOCF")) {
      strim <- c()
      for (j in seq(1, 8, by = 2)) {
        strim <- c(strim, paste0(all.resi[m, j], "(", all.resi[m, j+1], ")"))
      }
      names(strim) <- names(all.resi)[seq(1, 8, by = 2)]
      strim <- c(strim, unlist(all.resi[m, 9:12]))
      tbli[[paste0(ind, "ind", m)]] <- strim
    }
  }
  tbli <- as.data.frame(tbli)
  names(out.names) <- rownames(tbli)
  for (j in rownames(tbli)) {
    cat(paste(out.names[j], paste(unlist(tbli[j,]), collapse = " & "), sep = " & "), "\\\\", "\n")
  }
}


