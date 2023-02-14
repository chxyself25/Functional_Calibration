########### visualization: code run locally #####################
library(reshape2)
library(ggplot2)
library(fdapace)

file_code <- "FSH_BMI_INCOME_AGE"
fix_var <- c("BMI", "intercept", "income_high", "age")
res0 <- readRDS(paste0("./fcart_mr_", file_code, "_fix_us_Li.rds"))
bsres1 <- readRDS(paste0("./fcart_mr_", file_code, "_bs100_fix_us_Li.rds"))
bsres2 <- readRDS(paste0("./fcart_mr_", file_code, "_bs200_fix_us_Li.rds"))
names(bsres2[[1]])[3:5] <- c("intercept", "age", "income_high")
bsres0 <- lapply(1:length(bsres1), function(v) {rbind(bsres1[[v]], bsres2[[v]])})
names(bsres0) <- names(bsres1)
#res0 <- readRDS("./RealData/DataResults/fcart_mr_GLUCOSE_HRM_fix_sd_us.rds")
#bsres0 <- readRDS("./RealData/DataResults/fcart_mr_GLUCOSE_HRM_bs200_fix_sd_us.rds")

v <- "TRIGRES"
res <- res0[[v]]; bsres <- bsres0[[v]]
#names(bsres)[3:5] <- c("intercept", "age", "income_high")
plt <- NULL
vars <- sort(c(v, fix_var))
for (m in c("FCART", "LOCF")) {
  subres <- subset(res, method == m)
  #subres$race1 <- subres$race1-subres$race4
  #subres$race2_3 <- subres$race2_3-subres$race4
  #subres$race5 <- subres$race5-subres$race4
  #time_rescale <- diff(range(subres$time))
  #min_time <- min(subres$time)
  grid <- subres$time
  df <- melt(data.frame(subres[, vars], year = grid, method = m), 
             id.vars = c('year', 'method'), value.name = "estimate")
  df <- df[order(df$variable, df$year), ]
  #df <- data.frame(estimate = res$value, variable = res$variable, time = (grid-min_time)/time_rescale, method = m)
  bts <- lapply(1:200, function(b) {
    x <- subset(bsres, method == m & bts == b)
    y <- matrix(NA, nrow = length(grid), ncol = length(vars))
    indx <- which(grid < max(x$time) & grid > min(x$time))
    y[indx,] <- ConvertSupport(fromGrid = x$time, toGrid = grid[indx], phi = x[, vars])
    y
  })
  bsd = c()
  for (i in 1:length(vars)) {
    bsd = c(bsd, apply(sapply(bts, function(x) {x[,i]}), 1, sd, na.rm = TRUE))
  }
  df$upper <- df$estimate + 1.96*bsd
  df$lower <- df$estimate - 1.96*bsd
  #high_low <- df$estimate[df$variable == "income_levelhigh"] - df$estimate[df$variable == "income_levellow"]
  #h_l_sd <- sqrt(apply(sapply(bts, function(x) {x[,which(vars == "income_levelhigh")]}), 1, var, na.rm = TRUE) + 
  #  apply(sapply(bts, function(x) {x[,which(vars == "income_levellow")]}), 1, var, na.rm = TRUE))
  #df <- rbind(df, data.frame(year = grid, method = m, variable = "high-low", estimate = high_low, 
  #                           upper = high_low + 1.96*h_l_sd, lower = high_low - 1.96*h_l_sd))
  plt <- rbind(plt, df)
}
rg <- range(plt$year[plt$method == "FCART"])
ggplot(plt) + geom_line(aes(year, estimate), color = "black", size = 0.5) +
  geom_line(aes(year, upper), color = "blue", size = 0.5, linetype = 2) +
  geom_line(aes(year, lower), color = "blue", size = 0.5, linetype = 2) + 
  xlab("time (rescaled)") + ylab("slope") + facet_grid(variable~method, scales = "free_y") + 
  scale_x_continuous(breaks = seq(0, rg[2], by = 1), limits = rg) +
  theme_bw() + theme(text = element_text(size = 13))
#ggsave(file = paste0("./RealData/twoxt_FCART_LOCF_", v, "_", fix_var, "_slope_bts_lnsd.pdf"), width = 7, height = 6)

for (v in unique(plt$variable)) {
  rg <- range(plt$year[plt$method == "FCART"])
  yrg <- c(min(plt$lower[plt$variable == v & plt$year <= rg[2]], na.rm= TRUE), 
           max(plt$upper[plt$variable == v & plt$year <= rg[2]], na.rm = TRUE))
  # for (m in unique(plt$method)) {
  m = "FCART"
  pltmv <- subset(plt, method == m & variable == v) 
  v_name <- ifelse(v == "income_high", "HIGH INCOME", toupper(v))
  ggplot(pltmv) + geom_line(aes(year, estimate), color = "black", size = 0.5) + 
    geom_line(aes(year, upper), color = "blue", size = 0.5, linetype = 2) +
    geom_line(aes(year, lower), color = "blue", size = 0.5, linetype = 2) + 
    ylab(v_name) + ylim(yrg) + 
    scale_x_continuous(breaks = seq(0, rg[2], by = 1), limits = rg) + 
    theme_bw() + theme(text = element_text(size = 13))
  ggsave(file = paste0("./twoxt_z_", m, "_", v, "_slope_bts_income_age.pdf"), width = 4, height = 3)
}




