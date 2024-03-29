######### functional setting ###############
# generate X by xi
X_v <- function(tmat, xi_z, v) {
  # different variables share common eigen functions
  psi_1 <- psi(tmat, k = 1)
  psi_2 <- psi(tmat, k = 2)
  psi_3 <- psi(tmat, k = 3)
  n <- nrow(xi_z)
  m <- ncol(tmat)
  mu <- mean_fun(tmat, v)
  if (n != nrow(tmat)) {
    cat("The xi and t do not match!\n")
  }
  res <- mu + sweep(psi_1, 1, xi_z[,1], "*") + sweep(psi_2, 1, xi_z[,2], "*") + sweep(psi_3, 1, xi_z[,3], "*")
  c(t(res))
}

# define eigen functions: k = 1,2,3
psi <- function(t, k) {
  (1/sqrt(5))*sin(pi*k*t/10)
}

# mean function
mean_fun <- function(tt, v) {
  tt+v*sin(tt)
}

######### Error structure setting (to be extended to multivariate case) ################
evar_fun <- function(tt, ind = FALSE) {
  if (ind) {
    covs <- diag(1.5, nrow = length(tt))
  } else {
    covs <- outer(tt, tt, function(t1,t2) {
      2^(-abs(t1-t2)/5)
    })}
  return(covs)
}


multi_async_sim1 <- function(n = 200, m = 5, dx = 2, dz = 2, cor.pcs = c(0.8, 0.6, 0.4), sd.pcs = c(2, sqrt(2), 1), bx, bz) {
  ## xi covariance matrix
  XX <- matrix(0, ncol = 3*dx, nrow = 3*dx)
  for (i in 1:(3*dx-1)) {
    vars <- floor((i-1)/3) + 1
    pcs <- (i-1) %% 3 + 1
    rep.row <- rep(0, 3-pcs)
    if (vars < dx) {
      rep.row <- append(c(sapply((vars+1):dx, function(x) {c(0, 0, ((cor.pcs[pcs])^(x-vars))*(sd.pcs[pcs]^2))})),
                        rep.row)
    }
    XX[i, (i+1):(3*dx)] <- rep.row
  }
  XX <- XX + t(XX)
  diag(XX) <- rep(sd.pcs^2, dx)
  ## generate FPC
  xi_mat <- mvrnorm(n, mu = rep(0, 3*dx), Sigma = XX)
  ## generate Z: assume independence, pre-defined variance
  z_mat <- mvrnorm(n, mu = rep(0, dz), Sigma = diag(c(0.8, 1.2)))
  # generate t(response) observations
  t_obs <- t(sapply(1:n, function(x) {sort(runif(m, min = 0, max = 10))}))
  # generate s(covariate) observations
  s_obs <- list()
  for (v in 1:dx) {
    s_obs[[v]] <- t(sapply(1:n, function(x) {sort(runif(m, min = 0, max = 10))}))
  }
  # true X matrix is on T domain for all vairables
  X_t <- sapply(1:dx, function(x) {X_v(t_obs, xi_mat[, (3*(x-1)+1):(3*x)], x)})
  # true X matrix on S domain for all variables
  X_s <- sapply(1:dx, function(x) {X_v(s_obs[[x]], xi_mat[, (3*(x-1)+1):(3*x)], x)})
  # observed W on S domain for all variables, measurement error is iid with mean zero and variance sde for all varaibles
  W_s <- X_s + sapply(1:dx, function(x) {rnorm(m*n, mean = 0, sd = 1)})
  # time-invariant variables
  Z_ti <- apply(z_mat, 2, function(x) {rep(x, each = m)})
  # observed Y: generated by true underlying X_t and Z_t
  # error terms 
  e <- apply(t_obs, 1, function(tt) {
    mvrnorm(n=1, mu = rep(0,m), Sigma = evar_fun(tt))
  })
  Y_t <- cbind(X_t, Z_ti)%*%c(bx, bz) + c(e)
  
  ### formalize a data frame with observation time points saved
  time_df <- data.frame(ID = rep(1:n, each = m), time_t = c(t(t_obs)))
  for (v in 1:dx) {
    time_df[[paste0("time_s_", v)]] <- c(t(s_obs[[v]]))
  }
  
  return(list('Y_t' = Y_t, 'W_s' = W_s, 'X_t' = X_t, 'Z' = Z_ti, 
              'time_grid' = time_df))
}


multi_async_sim2 <- function(n = 200, m = 5, dx = 2, dz = 2, sd.pcs = sqrt(c(4,3,2,1)), bx, bz) {
  # generate t(response) observations
  t_obs <- t(sapply(1:n, function(x) {sort(runif(m, min = 0, max = 10))}))
  # generate s(covariate) observations
  s_obs <- list()
  for (v in 1:dx) {
    s_obs[[v]] <- t(sapply(1:n, function(x) {sort(runif(m, min = 0, max = 10))}))
  }
  # generate principle component scores
  xi_mat <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(sd.pcs^2))
  ## generate Z: assume independence
  z_mat <- mvrnorm(n, mu = rep(0, dz), Sigma = diag(c(0.8, 1.2)))
  # eigenfunctions construction
  f.basis <- create.fourier.basis(rangeval = c(0, 10*dx), nbasis = 4)
  # observed X: on s domain and t domain
  X_s <- NULL
  X_t <- NULL
  for (v in 1:dx) {
    meanvt <- mean_fun(t_obs, v)
    meanvs <- mean_fun(s_obs[[v]], v)
    # random sign: -1 or 1
    #signv <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
    X_s <- cbind(X_s, unlist(lapply(1:n, function(x) {
      meanvs[x,] + eval.basis(s_obs[[v]][x,]+(v-1)*10, f.basis)[,2:5] %*% xi_mat[x,]})))
    X_t <- cbind(X_t, unlist(lapply(1:n, function(x) {
      meanvt[x,] + eval.basis(t_obs[x,]+(v-1)*10, f.basis)[,2:5] %*% xi_mat[x,]})))
  }
  # observed W on s domain, adding measurement error
  W_s <- X_s + sapply(1:dx, function(x) {rnorm(m*n, mean = 0, sd = 1)})
  # time-invariant variables
  Z_ti <- apply(z_mat, 2, function(x) {rep(x, each = m)})
  # observed Y at time t
  # error terms 
  e <- apply(t_obs, 1, function(tt) {
    mvrnorm(n=1, mu = rep(0,m), Sigma = evar_fun(tt))
  })
  Y_t <- cbind(X_t, Z_ti)%*%c(bx, bz) + c(e)
  
  ### formalize a data frame with observation time points saved
  time_df <- data.frame(ID = rep(1:n, each = m), time_t = c(t(t_obs)))
  for (v in 1:dx) {
    time_df[[paste0("time_s_", v)]] <- c(t(s_obs[[v]]))
  }
  
  return(list('Y_t' = Y_t, 'W_s' = W_s, 'X_t' = X_t, 'Z' = Z_ti, 
              'time_grid' = time_df))
}

