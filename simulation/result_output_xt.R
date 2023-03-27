############################################
######### summarize results ################
############################################
library(dplyr); library(ggplot2)
library(fdapace)
library(reshape2)
beta_fun <- function(tt, intercept = FALSE, setting = 1) {
  if (intercept) {
    if (setting == 1) {
      betat = 0.2*tt + 0.5
    } else {
      betat <- tt^(1/2)
    }
  } else {
    if (setting == 1) {
      betat <- sin(pi*tt/10)
    } else {
      betat <- sin(pi*tt/5) 
    }
  }
  return(betat)
}
aeint <- function(beta, b, time, truncate = TRUE) {
  s_idx <- which(!is.na(beta))
  if (truncate) {
    s_idx <- intersect(s_idx, which(time < 9.75 & time > 0.25))
  }
  res <- trapzRcpp(time[s_idx], abs((beta-b)[s_idx]))
  
  return(res)
}
seint <- function(beta, b, time, truncate = TRUE) {
  s_idx <- which(!is.na(beta))
  if (truncate) {
    s_idx <- intersect(s_idx, which(time < 9.75 & time > 0.25))
  }
  res <- trapzRcpp(time[s_idx], ((beta-b)[s_idx])^2)
  return(res)
}
interp_NA <- function(f, t, y) {
  res = rep(NA, length(t))
  if (all(is.na(y))) {
    return(res)
  }
  f <- f[!is.na(y)]
  y <- y[!is.na(y)]
  idx <- which(t <= max(f) & t >= min(f))
  if (length(idx) == 0) {
    return(res)
  }
  res[idx] <- ConvertSupport(f, t[idx], y)
  return(res)
}

m = "m15"
for (set in 1:2) {
  res <- readRDS(paste(c(paste0("./onext", set, "_sim_nind"), m, "200.rds"), collapse = "_"))
  if (!is.null(m) && m == "m15" & set == 2) {
    id_seq <- unlist(lapply(1:10, function(i) {sapply(1:20, function(iter) {iter + 20*(i-1)})}))
    res$id <- rep(id_seq, each = 4*60) 
  }
  res$method <- gsub("FM", "FCAR", toupper(res$method))
  res$method <- gsub("\\<F\\>", "FVCM", res$method)
  res$method <- gsub("\\<O\\>", "Oracle", res$method)
  res$method <- gsub("CZF", "KW", res$method)
  #res <- subset(res, method != "CZF")
  #czf <- readRDS("./OneX/Simulations/onext1_sim_nind_czf_200.rds")
  #res <- bind_rows(res, czf)
  res$b0 = beta_fun(res$time, intercept = TRUE, set)
  res$b1 = beta_fun(res$time, intercept = FALSE, set)
  #res <- subset(res, !grepl("int", metric))
  #res$metric <- as.character(res$metric); res$method <- as.character(res$method)
  res2 <- res %>%
    group_by(method, id) %>%
    summarize(
      MADE = (1/20)*(aeint(beta0, b0, time)/diff(beta_fun(c(0,10),intercept = TRUE,set))+
                       aeint(beta1, b1, time)/diff(range(beta_fun(seq(0,10,by=0.01),FALSE,set)))),
      WASE = (1/20)*(seint(beta0, b0, time)/diff(beta_fun(c(0,10),intercept = TRUE,set))^2+
                       seint(beta1, b1, time)/diff(range(beta_fun(seq(0,10,by=0.01),FALSE,set)))^2)
    ) %>% as.data.frame
  df <- melt(res2, id = c("method", "id"), variable.name = "metric", value.name = "value") %>% 
    group_by(method, metric) %>%
    summarise(mean = round(mean(value), 3), sd= round(sd(value), 3),
              med = round(median(value), 3), 
              p25 = unname(round(quantile(value, 0.25), 3)), 
              p75 = unname(round(quantile(value, 0.75), 3))) %>% as.data.frame
  df$metric <- as.character(df$metric)
  header1 = ifelse(set == 1, "Setting I", "Setting \\II")
  for (c in c("MADE", "WASE")) {
    header2 <- c
    for (md in c("FCAR", "FVCM", "KW", "Oracle")) {
      out_str <- paste(c(header1, md, header2), collapse = " & ")
      df_sub <- subset(df, method == md & metric == c)
      mean_sd <- paste0(df_sub$mean[1], "(", df_sub$sd[1], ")")
      out_str <- paste(c(out_str, mean_sd, df_sub[, c("med", "p25", "p75")]), collapse = " & ")
      cat(out_str, "\\\\", "\n")
      header1 = ""; header2 = ""
    } 
    cat("\\hline", "\n")
  }
}
# df$metric <- as.character(df$metric)
# for (m in unique(df$metric)) {
#   dfm <- subset(df, metric == m)
#   str <- c()
#   for (j in 1:5) {
#     strj <- c(dfm[j, 1:2], paste0(dfm[j,3], "(", dfm[j,4], ")"), unlist(dfm[j,5:7]))
#     str <- append(str, paste(strj, collapse = " & "))
#   }
#   print(paste(str, collapse = " \\ "))
# }


############################################
######## summarize results in plots ########
############################################
percent <- 100
grid <- seq(0+0.1*(100-percent)*0.5, 10-0.1*(100-percent)*0.5, by = 0.1)
m = NULL
for (set in 1:2) {
  res <- readRDS(paste(c(paste0("./onext", set, "_sim_nind"), m, "200.rds"), collapse = "_"))
  if (!is.null(m) && m == "m15" & set == 2) {
    id_seq <- unlist(lapply(1:10, function(i) {sapply(1:20, function(iter) {iter + 20*(i-1)})}))
    res$id <- rep(id_seq, each = 4*60) 
  }
  res$method <- gsub("FM", "FCAR", toupper(res$method))
  res$method <- gsub("\\<F\\>", "FVCM", res$method)
  res$method <- gsub("\\<O\\>", "Oracle", res$method)
  res$method <- gsub("CZF", "KW", res$method)
  
  res30 <- res %>%
    group_by(method, id) %>%
    summarize(beta0 = interp_NA(time, grid, beta0), time = grid) %>%
    ungroup() %>%
    group_by(method, time) %>%
    summarize(avg = median(beta0, na.rm = TRUE), upper = quantile(beta0, 0.975, na.rm=TRUE), 
              lower = quantile(beta0, 0.025, na.rm=TRUE)) %>% 
    mutate(true0 = beta_fun(time, intercept = TRUE, set)) %>% as.data.frame
  ggplot(res30) + geom_line(aes(time, avg), color = "black", size = 0.5) +
    geom_line(aes(time, true0), color = "red", size = 0.5) +
    geom_line(aes(time, upper), color = "blue", size = 0.5, linetype = 2) +
    geom_line(aes(time, lower), color = "blue", size = 0.5, linetype = 2) +
    facet_wrap(.~method, nrow = 2, ncol = 3, scale = "free_y") + theme_bw() + ylab("b0")
  #ggsave(paste0("./OneX/all5_beta0_sim_xt", set, "_bands", bw, ".pdf"), width = 6, height = 5)
  rng <- c(min(res30$lower[res30$method %in% c("FCAR", "Oracle")], na.rm = TRUE), max(res30$upper[res30$method %in% c("FCAR", "Oracle")], na.rm = TRUE))
  for (md in unique(res30$method)) {
    res30m <- subset(res30, method == md)
    plt <- ggplot(res30m) + geom_line(aes(time, avg), color = "black", size = 0.8) + 
      geom_line(aes(time, true0), color = "red", size = 0.8) + 
      geom_line(aes(time, upper), color = "blue", size = 0.8, linetype = 2) + 
      geom_line(aes(time, lower), color = "blue", size = 0.8, linetype = 2) + 
      ylab("intercept") + theme_bw() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
    if (md %in% c("FCAR", "Oracle")) {
      plt <- plt + ylim(rng)
    }
    ggsave(plt, file = paste(c(paste0("./", md,"_beta0_sim_xt", set), m, "bands", paste0(percent, ".pdf")), collapse = "_"))
  }
  
  res31 <- res %>%
    group_by(method, id) %>%
    summarize(beta1 = interp_NA(time, grid, beta1), time = grid) %>%
    ungroup() %>%
    group_by(method, time) %>%
    summarize(avg = median(beta1, na.rm = TRUE), upper = quantile(beta1, 0.975, na.rm=TRUE), 
              lower = quantile(beta1, 0.025, na.rm=TRUE)) %>% 
    mutate(true1 = beta_fun(time, intercept = FALSE, set)) %>% as.data.frame
  ggplot(res31) + geom_line(aes(time, avg), color = "black", size = 0.5) +
    geom_line(aes(time, true1), color = "red", size = 0.5) +
    geom_line(aes(time, upper), color = "blue", size = 0.5, linetype = 2) +
    geom_line(aes(time, lower), color = "blue", size = 0.5, linetype = 2) +
    facet_wrap(.~method, nrow = 2, ncol = 3, scale = "free_y") + theme_bw() + ylab("b1")
  rng <- c(min(res31$lower[res31$method %in% c("FCAR", "Oracle")], na.rm = TRUE), max(res31$upper[res31$method %in% c("FCAR", "Oracle")], na.rm = TRUE))
  for (md in unique(res31$method)) {
    res31m <- subset(res31, method == md)
    plt <- ggplot(res31m) + geom_line(aes(time, avg), color = "black", linetype = "longdash", size = 0.8) + 
      geom_line(aes(time, true1), color = "black", size = 0.8) + 
      geom_line(aes(time, upper), color = "blue", size = 0.8, linetype = 2) + 
      geom_line(aes(time, lower), color = "blue", size = 0.8, linetype = 2) + 
      ylab("slope") + theme_bw() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
    if (md %in% c("FCAR", "Oracle")) {
      plt <- plt + ylim(rng)
    }
    #ggsave(paste0("./", md,"_beta1_sim_xt", set, "_bands", bw, ".pdf"))
    ggsave(plt, file = paste(c(paste0("./", md,"_beta1_sim_xt", set), m, "bands", paste0(percent, ".pdf")), collapse = "_"))
  }
}


##### computation time comparison ###########
res_all <- NULL
m = "noerror"
for (set in 1:2) {
  res <- readRDS(paste(c(paste0("./onext", set, "_sim_nind"), m, "200.rds"), collapse = "_"))
  if (m == "m15" & set == 2) {
    id_seq <- unlist(lapply(1:10, function(i) {sapply(1:20, function(iter) {iter + 20*(i-1)})}))
    res$id <- rep(id_seq, each = 4*60) 
  }
  res_all <- rbind(res_all, res %>% mutate(setting = set))
}
df <- res_all %>% group_by(method) %>% 
  summarize(mean_time = mean(run_time), ct = n()) %>% as.data.frame
df
