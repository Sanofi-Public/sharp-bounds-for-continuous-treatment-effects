# Functions used to simulate the dataset


simulated_dataset <- function(
  # Sample size
  n=1e3,
  # Number of confounders
  p_X=5,
  p_U=5,
  # Correlation between confounders
  rho_X=0.1,
  rho_U=0.1,
  # Proportion of the maximum correlation between X and U to take
  corr_XU_prop=0.6,
  # Effects on treatment
  beta_X=rep(0.2, p_X),
  beta_U=rep(0.2, p_U),
  # Effects on response
  gamma_X=rep(0.25, p_X),
  gamma_U=rep(0.2, p_U), 
  zeta=-1,
  # Standard deviations
  sd_eps_T=NULL,
  sd_eps_Y=NULL) {
  
  stopifnot(0 <= rho_X, rho_X <= 1)
  stopifnot(0 <= rho_U, rho_U <= 1)
  stopifnot(0 <= corr_XU_prop, corr_XU_prop <= 1)
  
  # Covariance matrices
  Sigma_X <- matrix(diag(rep(1, p_X)), ncol=p_X)
  Sigma_X[abs(row(Sigma_X)-col(Sigma_X)) == 1] <- rho_X
  
  Sigma_U <- matrix(diag(rep(1, p_U)), ncol=p_U)
  Sigma_U[abs(row(Sigma_U)-col(Sigma_U)) == 1] <- rho_U
  
  # Maximum value for the correlation between X and U to have a diagonal dominant matrix (invertible)
  corr_XU_max <- (1 - rho_X) / p_U
  corr_XU <- corr_XU_max * corr_XU_prop
  cov_XU <- matrix(corr_XU, nrow=p_X, ncol=p_U)
  Sigma <- rbind(cbind(Sigma_X, cov_XU),
                 cbind(t(cov_XU), Sigma_U))
  
  # Simulation of (X, U)
  XU <- mult.normal(n, mu=rep(0, p_X+p_U), Sigma)
  X <- XU[, 1:p_X]
  U <- XU[, (p_X+1):(p_X+p_U)]
  
  # Simulation of T
  mu_T_XU <- XU %*% c(beta_X, beta_U)
  if (is.null(sd_eps_Y) & is.null(sd_eps_T)) {
    sd_eps_Y <- sd_eps_T <- sd(mu_T_XU) / 4  # Here, 4 is arbitrary
  }
  
  epsilon <- rnorm(n, 0, sd_eps_T)
  T_XU <- mu_T_XU + epsilon  # Simulated treatment

  Y_XU <- T_XU + zeta * (X %*% gamma_X) * exp(- T_XU * (X %*% gamma_X)) - (U %*% gamma_U) * (X %*% gamma_X) + rnorm(n, 0, sd_eps_Y)
  
  # Removing gross outliers
  hat_values <- hatvalues(lm(rnorm(n, 0, 1) ~ T_XU + X + Y_XU))
  outliers <- which(hat_values > quantile(hat_values, 0.9))
  n_effective <- length(T_XU[-outliers])  # Effective sample size after removing outliers
  
  data <- as.data.frame(X[-outliers, ])
  data$t <- T_XU[-outliers]
  data$Y <- Y_XU[-outliers]
  
  simu <- NULL
  simu$data <- data
  simu$mu_T_XU <- mu_T_XU
  simu$p_X <- p_X
  simu$p_U <- p_U
  simu$rho_X <- rho_X
  simu$rho_U <- rho_U
  simu$beta_X <- beta_X
  simu$beta_U <- beta_U 
  simu$gamma_X <- gamma_X 
  simu$gamma_U <- gamma_U 
  simu$zeta <- zeta 
  simu$sd_eps_Y <- sd_eps_Y
  simu$sd_eps_T <- sd_eps_T
  simu$Sigma_X <- Sigma_X
  simu$Sigma_U <- Sigma_U
  simu$cov_XU <- cov_XU
  simu$outliers <- outliers
  simu$n_effective <- n_effective
  
  return(simu)
}


mu_U_x <- function(cov_XU, Sigma_X, X) {
  return(t(t(cov_XU) %*% solve(Sigma_X) %*% t(X)))
}


sigma_U_x <- function(cov_XU, Sigma_X, Sigma_U) {
  return(Sigma_U - t(cov_XU) %*% solve(Sigma_X) %*% cov_XU)
}


sigma_X_u <- function(cov_XU, Sigma_X, Sigma_U) {
  return(Sigma_X - cov_XU %*% solve(Sigma_U) %*% t(cov_XU))
}


mu_T_x <- function(X, beta_X, beta_U, cov_XU, Sigma_X) {
  mu_cond <- X %*% beta_X + mu_U_x(cov_XU, Sigma_X, X) %*% beta_U
  return(t(mu_cond))
}


sigma_T_x <- function(sigma_eps, beta_U, cov_XU, Sigma_X, Sigma_U) {
  var_cond <- sigma_eps^2 + t(beta_U) %*% sigma_U_x(cov_XU, Sigma_X, Sigma_U) %*% beta_U
  return(var_cond)
}


part_explained_var <- function(sigma_eps, beta_X, cov_XU, Sigma_X, Sigma_U) {
  sigma_T_u <- t(beta_U) %*% sigma_X_u(cov_XU, Sigma_X, Sigma_U) %*% beta_U
  sigma_T <- sigma_eps^2 + t(beta_X) %*% Sigma_X %*% beta_X + t(beta_U) %*% Sigma_U %*% beta_U + t(beta_X) %*% cov_XU %*% beta_U
  return(sigma_T_u / sigma_T)
}


mu_Y_t_XU <- function(X, U, gamma_X, gamma_U, t, zeta) {
  mean_cond <- t - zeta * (X %*% gamma_X) * exp(- t * (X %*% gamma_X)) - (U %*% gamma_U) * (X %*% gamma_X)
  return(mean_cond)
}


capo_t_X <- function(X, gamma_X, gamma_U, t, zeta, cov_XU, Sigma_X) {
  X_gamma_X <- X %*% gamma_X
  mu_Y_t_X <- t + zeta * X_gamma_X * exp(- t * X_gamma_X) - (mu_U_x %*% gamma_U) * X_gamma_X
  return(mu_Y_t_X)
}


apo_t <- function(t, zeta, gamma_X, gamma_U, cov_XU, Sigma_X) {
  sig_x <- c(t(gamma_X) %*% Sigma_X %*% gamma_X)
  mu_Y_t <- t * ( 1 - zeta * sig_x * exp((t^2 / 2) * sig_x)) - c(t(gamma_U) %*% t(cov_XU) %*% gamma_X)
  return(mu_Y_t)
}


mult.normal <- function(n, mu, Sigma) {
  p <- length(mu)
  eigen.val.vec <- eigen(Sigma, symmetric=TRUE)
  X <- drop(mu) + eigen.val.vec$vectors %*% diag(sqrt(pmax(eigen.val.vec$values, 0)), p) %*% t(matrix(rnorm(p*n), n))
  mu.names <- names(mu)
  Sigma.names <- dimnames(Sigma)
  
  if (!is.null(Sigma.names) && is.null(mu.names)) {
    mu.names <- Sigma.names[[1]]
  }
  
  dimnames(X) <- list(mu.names, NULL)
  
  return(t(X))
}

