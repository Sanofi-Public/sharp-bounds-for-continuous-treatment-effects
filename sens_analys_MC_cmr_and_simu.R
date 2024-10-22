# In this file, we perform the sensitivity analyses on each Monte-Carlo simulated sample and on real data

# Import libraries
library(foreach)
library(torch)
library(latex2exp)  # For LaTeX expressions
library(tictoc)

# Import R files
source("./utils.R")

# Get provided arguments
args <- commandArgs()
# 1: parallel.computation.jesson (boolean)
# 2: n.MC (1:infty)
# 3: data.name ("simul" or "cmr")
# 4: gamma_est, for cmr data (real number from 1 to 50) 
# 5: xi.method (NULL or "neural_network")
# 6: compute.jesson (boolean)
# 7: doses.length (2:infty)
# 8: version (1:infty)
parallel.computation.jesson.arg <- as.logical(args[1])
n.MC.arg <- as.numeric(args[2])
data.name.arg <- args[3]
gamma_est.arg <- as.numeric(args[4])
if (args[5] == "NULL") {
  xi.method.arg <- NULL
} else {
  xi.method.arg <- args[5]
}
compute.jesson.arg <- as.logical(args[6])
doses.length.arg <- as.numeric(args[7])
version.arg <- as.numeric(args[8])

# Whether to show or not progress bars
verbose <- FALSE

if (cuda_is_available()) {
  nb.cuda.device <- cuda_device_count()
  message(paste("Number of CUDA devices:", nb.cuda.device))
  device <- torch_device("cuda")
} else {
  message("CUDA unavailable")
  device <- torch_device("cpu")
}

seed <- 1
set.seed(seed)
torch_manual_seed(seed)


parallel.computation.jesson <- parallel.computation.jesson.arg  # FALSE  # TRUE

if (parallel.computation.jesson) {
  nb.cores <- parallel::detectCores() - 1  # Get number of cores and subtract 1 or 2 (DO NOT USE ALL THE CORES)
  cl <- parallel::makeCluster(nb.cores)  # Create the cluster
  invisible(parallel::clusterExport(cl, c("grid.search.optimizer", "CAPO.estim", "heaviside.fun")))
  doSNOW::registerDoSNOW(cl)  # Register the cluster
}


# Number of Monte-Carlo samples
n.MC <- n.MC.arg  # 1  # 20  # 50

doses <- NULL

# Lists to store the results on each Monte-Carlo sample
qb.mc.bounds <- list()
jesson.mc.bounds <- list()

# Lists to store the execution times on each Monte-Carlo sample
qb.exec.times <- rep(NA, n.MC)
jesson.exec.times <- rep(NA, n.MC)

# Gamma used for each Monte-Carlo sample. Finally, each Monte-Carlo sample will use the same sensitivity parameter Gamma
gammas <- rep(NA, n.MC)

for (mc.ind in 1:n.MC) {
  
  message(paste("Monte-Carlo sample:", mc.ind, "/", n.MC))
  
  data.name <- data.name.arg  # "simul"  # "cmr"
  
  if (data.name == "cmr") {
    
    # Folder name to retrieve the data
    data.folder.name <- "data/pm25/"
    
    # Preprocess the data
    preprocessed.data <- preprocess.pm2.5.cmr.data(data.folder.name)
    
    unnormalized.data <- preprocessed.data$unnormalized.data  # Get all unnormalized data
    all.data <- preprocessed.data$normalized.data  # Get all preprocessed data.frame
    n.all <- preprocessed.data$n.all  # Get number of individuals
    cov.names <- preprocessed.data$cov.names  # Get covariate names
    scaled.Y <- preprocessed.data$scaled.Y  # Get scaled outcome (CMR)
    scaled.t <- preprocessed.data$scaled.t  # Get scaled treatment (PM2.5)
    
    X <- subset(all.data, select=-c(Y, t))
    t <- all.data$t
    Y <- all.data$Y
    
    n <- n.effective <- nrow(X)
    
    gamma_est <- gamma_est.arg  # Value of the sensitivity parameter to use
    quant.order <- NULL
    
  } else if (data.name == "simul") {
    
    ### Simulation
    
    # Sample size
    n <- 1000
    
    # Number of confounders
    p_X <- 5
    p_U <- 3
    
    # Correlation between confounders
    rho_X <- 0.3
    rho_U <- 0.3
    corr_XU_prop <- 0.5
    beta_X <- rep(0.3, p_X)
    beta_U <- rep(0.2, p_U)
    
    gamma_X <- rep(0.2, p_X)
    gamma_U <- c(rep(0.4, floor(p_U/2)), rep(0.7, p_U-floor(p_U/2)))
    zeta <- -0.3
    # Observed Y
    sd_eps_T <- 0.5
    sd_eps_Y <- 0.3
    
    if (mc.ind == 1) {  # Estimate Gamma on the first dataset that we put then aside
      
      simu <- simulated_dataset(n=n, p_X=p_X, p_U=p_U,
                                rho_X=rho_X, rho_U=rho_U,
                                corr_XU_prop=corr_XU_prop,
                                beta_X=beta_X, beta_U=beta_U,
                                gamma_X=gamma_X, gamma_U=gamma_U,
                                zeta=zeta, sd_eps_T=sd_eps_T, sd_eps_Y=sd_eps_Y)
      
      n.effective <- nrow(simu$data)
      
      grid_t <- seq(min(simu$data$t), max(simu$data$t), 0.1)
      cov_XU <- simu$cov_XU
      Sigma_X <- simu$Sigma_X

      X <- subset(simu$data, select=-c(Y, t))
      t <- simu$data$t
      Y <- simu$data$Y
      
      # Estimation of the sensitivity parameter
      pdf_T_XU <- dnorm(t, mean=simu$mu_T_XU[-simu$outliers], sd=simu$sd_eps_T) + 1e-4
      pdf_T_X <- dnorm(t, mean=mu_T_x(as.matrix(X), simu$beta_X, simu$beta_U, simu$cov_XU, simu$Sigma_X), sd=sqrt(sigma_T_x(simu$sd_eps_T, simu$beta_U, simu$cov_XU, simu$Sigma_X, simu$Sigma_U))) + 1e-4
      
      quotient <- pdf_T_XU / pdf_T_X
      quant.order <- 0.99
      gamma_est <- unname(quantile(quotient, probs=c(quant.order)))
      
    }
    
    simu <- simulated_dataset(n=n, p_X=p_X, p_U=p_U,
                              rho_X=rho_X, rho_U=rho_U,
                              corr_XU_prop=corr_XU_prop,
                              beta_X=beta_X, beta_U=beta_U,
                              gamma_X=gamma_X, gamma_U=gamma_U,
                              zeta=zeta, sd_eps_T=sd_eps_T, sd_eps_Y=sd_eps_Y)
    
    n.effective <- nrow(simu$data)
    
    grid_t <- seq(min(simu$data$t), max(simu$data$t), 0.1)
    cov_XU <- simu$cov_XU
    Sigma_X <- simu$Sigma_X

    X <- subset(simu$data, select=-c(Y, t))
    t <- simu$data$t
    Y <- simu$data$Y
    
    apo_for_plot_fun <- function(t) apo_t(t=t, zeta=zeta,
                                          gamma_X=gamma_X, gamma_U=gamma_U,
                                          cov_XU=cov_XU, Sigma_X=Sigma_X)
    
  } else {
    stop("data.name must be 'cmr' or 'simul'")
  }
  
  
  ### Fine-tuning
  
  lr.space <- seq(from=0.0001, to=0.001, by=0.0001)
  K.space <- seq(from=3, to=30, by=1)
  dim.hidden.space <- c(8, 16, 32, 64)
  
  # Number of triplets of hyperparameters to test
  n.hyper <- 100
  
  rd.lr.ind <- sample(1:length(lr.space), n.hyper, replace=TRUE)
  rd.K.ind <- sample(1:length(K.space), n.hyper, replace=TRUE)
  rd.dim.hidden.ind <- sample(1:length(dim.hidden.space), n.hyper, replace=TRUE)
  
  lr.vec <- lr.space[rd.lr.ind]
  K.vec <- K.space[rd.K.ind]
  dim.hidden.vec <- dim.hidden.space[rd.dim.hidden.ind]
  
  # Proportion of data from dataset D in D1
  D1.prop <- 0.5
  
  # Parameters for hyperparameters fine-tuning
  n.random.splits <- 2
  max.iter <- 2000
  patience <- 20
  
  train.prop <- 0.8
  valid.prop <- 0.1
  
  # Parameters for our method
  # Range of windows to test
  bandwidths.length <- 40
  bandwidths <- seq(from=0.1, to=2.5, length.out=bandwidths.length)
  use.stabilization <- TRUE
  cond.quant.method <- "root_search"  # "quantile_forest"
  retrain.cond.quant <- FALSE  # TRUE
  xi.method <- xi.method.arg  #NULL  # "neural_network"  # "regression_forest"
  
  # Parameters for Jesson et al.'s method
  # num.trees <- 2000  # Number of trees used for regression_forest
  Y.sample.len <- 500  # Size of the sample of observed outcomes used in the grid search
  X.sample.len <- NULL  # Number of covariates on which to take the mean to compute the APO bounds
  
  compute.CI <- TRUE
  
  # Parameters for confidence intervals
  
  # Confidence interval level
  alpha <- 0.05
  alpha.low <- alpha / 2
  alpha.high <- 1 - alpha / 2
  
  # Number of bootstrap resamples
  B <- 100
  
  max.iter.gps <- 1000
  max.iter.resp <- 1000
  patience.gps <- 5
  patience.resp <- 5
  
  fine.tun.nn.params <- list(train.prop.gps=0.8, valid.prop.gps=0.1, n.random.splits.gps=2, max.iter.gps=max.iter.gps, patience.gps=patience.gps,
                             train.prop.resp=0.8, valid.prop.resp=0.1, n.random.splits.resp=2, max.iter.resp=max.iter.resp, patience.resp=patience.resp)
  nn.params <- list(max.iter.gps=max.iter.gps, patience.gps=patience.gps,
                    max.iter.resp=max.iter.resp, patience.resp=patience.resp)
  nn.init <- NULL
  
  compute.qb <- TRUE
  compute.jesson <- compute.jesson.arg  # FALSE  # TRUE
  
  # Sensitivity parameter
  gamma <- gamma_est
  gammas[mc.ind] <- gamma

  if (is.null(doses)) {
    
    # Exposures for which we want to compute the intervals
    doses.length <- doses.length.arg  # 5  # 15
    doses <- seq(from=as.numeric(quantile(t, probs=0.05)),
                 to=as.numeric(quantile(t, probs=0.95)),
                 length.out=doses.length)
    
  }
  
  
  # Start measuring execution time
  tic("Sensitivity analysis")
  
  #### PEI computation on D2 ####
  message("----- PEI estimation -----")
  
  qb.jesson.pei <- PEI.fun(X=X, Y=Y, t=t,
                           data.name=data.name,
                           dose=doses,
                           gamma=gamma,
                           bootstrap.ind=NULL,
                           stabilization=use.stabilization,
                           bandwidths=bandwidths,
                           B.param=NULL,
                           compute.QB=compute.qb,
                           cond.quant.method=cond.quant.method,
                           Q.predict=NULL,
                           xi.method=xi.method,
                           #retrain.cond.quant=retrain.cond.quant,
                           compute.Jesson=compute.jesson,
                           X.sample.len=NULL,
                           Y.sample.len=Y.sample.len,
                           nn.init=NULL,
                           xi.models=NULL,
                           D1.prop=D1.prop,
                           fine.tun.nn.params=fine.tun.nn.params,
                           nn.params=nn.params,
                           grid.K=K.vec,
                           grid.hid.dim=dim.hidden.vec,
                           grid.lr=lr.vec,
                           use.parallel.jesson=parallel.computation.jesson,
                           device=device,
                           verbose=verbose)
  
  Q.predict <- qb.jesson.pei$Q.predict
  nn.init <- qb.jesson.pei$nn.init  # This can be used after when d >= 2 not to search again the best hyperparameters
  xi.models <- qb.jesson.pei$xi.models
  jesson.unconf.vec <- qb.jesson.pei$jesson.unconf.per.dose
  exec.time.qb.PEI <- qb.jesson.pei$exec.time.qb
  exec.time.jesson.PEI <- qb.jesson.pei$exec.time.jesson
  
  # Prepare the neural network hyperparameters for the bootstrap samples
  nn.params <- list(max.iter.gps=max.iter.gps, patience.gps=patience.gps,
                    K.gps=qb.jesson.pei$optimal.params.gps$K.optim,
                    lr.gps=qb.jesson.pei$optimal.params.gps$lr.optim,
                    hid.dim.gps=qb.jesson.pei$optimal.params.gps$dim.hidden.optim,
                    max.iter.resp=max.iter.resp, patience.resp=patience.resp,
                    K.resp=qb.jesson.pei$optimal.params.resp$K.optim,
                    lr.resp=qb.jesson.pei$optimal.params.resp$lr.optim,
                    hid.dim.resp=qb.jesson.pei$optimal.params.resp$dim.hidden.optim)
  
  # Compute the CI
  qb.boot.lb.list <- list()
  qb.boot.ub.list <- list()
  
  jesson.boot.lb.mat <- matrix(NA, nrow=B, ncol=doses.length)
  jesson.boot.ub.mat <- matrix(NA, nrow=B, ncol=doses.length)
  
  exec.times.qb.vec <- rep(NA, B)
  exec.times.jesson.vec <- rep(NA, B)
  
  # Save the PEI for each bootstrap resample (I do not use jesson.boot.pei)
  jesson.boot.pei <- foreach(b=1:B) %do% {
    
    message(paste("##### Bootstrap resample:", b, "/", B, "#####"))
    
    # Non parametric bootstrap
    bootstrap.ind <- sample(1:n.effective, n.effective, replace=TRUE)
    # Get bootstrap covariates
    X.boot <- X[bootstrap.ind, ]
    t.boot <- t[bootstrap.ind]
    Y.boot <- Y[bootstrap.ind]
    
    qb.jesson.pei.boot <- PEI.fun(X=X.boot, Y=Y.boot, t=t.boot,
                                  data.name=data.name,
                                  dose=doses,
                                  gamma=gamma,
                                  bootstrap.ind=bootstrap.ind,
                                  stabilization=use.stabilization,
                                  bandwidths=bandwidths,
                                  B.param=NULL,
                                  compute.QB=compute.qb,
                                  cond.quant.method=cond.quant.method,
                                  xi.method=xi.method,
                                  #retrain.cond.quant=retrain.cond.quant,
                                  Q.predict=Q.predict[bootstrap.ind, ],  # Reorder the conditional quantiles
                                  compute.Jesson=compute.jesson,
                                  X.sample.len=NULL,
                                  Y.sample.len=Y.sample.len,
                                  nn.init=nn.init,
                                  xi.models=xi.models,
                                  D1.prop=D1.prop,
                                  fine.tun.nn.params=NULL,
                                  nn.params=nn.params,
                                  grid.K=NULL,
                                  grid.hid.dim=NULL,
                                  grid.lr=NULL,
                                  use.parallel.jesson=parallel.computation.jesson,
                                  device=device,
                                  verbose=verbose)
    
    qb.boot.lb.list[[b]] <- qb.jesson.pei.boot$qb.PEI.per.window.and.dose.lb
    qb.boot.ub.list[[b]] <- qb.jesson.pei.boot$qb.PEI.per.window.and.dose.ub
    
    jesson.boot.lb.mat[b, ] <- qb.jesson.pei.boot$jesson.PEI.per.dose.lb
    jesson.boot.ub.mat[b, ] <- qb.jesson.pei.boot$jesson.PEI.per.dose.ub
    
    exec.times.qb.vec[b] <- qb.jesson.pei.boot$exec.time.qb
    exec.times.jesson.vec[b] <- qb.jesson.pei.boot$exec.time.jesson
  }
  
  start.time.qb1 <- Sys.time()
  
  # Compute the MSE for each dose and each bandwidth
  qb.mse.lb <- matrix(0, nrow=bandwidths.length, ncol=doses.length)
  qb.mse.ub <- matrix(0, nrow=bandwidths.length, ncol=doses.length)
  
  for (b in 1:B) {
    
    new.term.lb <- (qb.boot.lb.list[[b]] - qb.jesson.pei$qb.PEI.per.window.and.dose.lb)**2
    new.term.ub <- (qb.boot.ub.list[[b]] - qb.jesson.pei$qb.PEI.per.window.and.dose.ub)**2
    
    qb.mse.lb <- qb.mse.lb + new.term.lb
    qb.mse.ub <- qb.mse.ub + new.term.ub
    
  }
  
  # If NA value, put Inf instead
  qb.mse.lb[is.na(qb.mse.lb)] <- Inf
  qb.mse.ub[is.na(qb.mse.ub)] <- Inf
  
  # This step is useless if we only want to find the index of the minimum
  qb.mse.lb <- qb.mse.lb / B
  qb.mse.ub <- qb.mse.ub / B
  
  # Get index of optimal windows h- and h+ for our method
  h.lb.ind <- apply(FUN=which.min, X=qb.mse.lb, MARGIN=2)
  h.ub.ind <- apply(FUN=which.min, X=qb.mse.ub, MARGIN=2)
  
  if (doses.length == 1) {  # To be tested
    
    # Store the optimal PEIs for continuous QB
    qb.PEI.lb.vec <- qb.jesson.pei$qb.PEI.per.window.and.dose.lb[h.lb.ind]
    qb.PEI.lb.vec <- qb.PEI.lb.vec[1, 1]
    qb.PEI.ub.vec <- qb.jesson.pei$qb.PEI.per.window.and.dose.ub[h.ub.ind]
    qb.PEI.ub.vec <- qb.PEI.ub.vec[1, 1]
    
    qb.unconf.lb.vec <- qb.jesson.pei$qb.unconf.per.window.and.dose.lb[h.lb.ind, ]
    qb.unconf.lb.vec <- qb.unconf.lb.vec[1, 1]
    qb.unconf.ub.vec <- qb.jesson.pei$qb.unconf.per.window.and.dose.ub[h.ub.ind, ]
    qb.unconf.ub.vec <- qb.unconf.ub.vec[1, 1]
    qb.unconf.vec <- (qb.unconf.lb.vec + qb.unconf.ub.vec) / 2  # An average because the chosen h is different between lower and upper bounds
    
  } else {
    
    # Store the optimal PEIs for continuous QB
    qb.PEI.lb.vec <- diag(qb.jesson.pei$qb.PEI.per.window.and.dose.lb[h.lb.ind, ])
    qb.PEI.ub.vec <- diag(qb.jesson.pei$qb.PEI.per.window.and.dose.ub[h.ub.ind, ])
    
    qb.unconf.lb.vec <- diag(qb.jesson.pei$qb.unconf.per.window.and.dose.lb[h.lb.ind, ])
    qb.unconf.ub.vec <- diag(qb.jesson.pei$qb.unconf.per.window.and.dose.ub[h.ub.ind, ])
    qb.unconf.vec <- (qb.unconf.lb.vec + qb.unconf.ub.vec) / 2  # An average because the chosen h is different between lower and upper bounds
    
  }
  
  # Compute the CIs
  # For continuous QB
  qb.lb.ub.optimal.mat <- foreach(b=1:B, .combine="rbind") %do% {
    c(diag(qb.boot.lb.list[[b]][h.lb.ind, ]), diag(qb.boot.ub.list[[b]][h.ub.ind, ]))
  }
  
  qb.CI.lb.vec <- apply(FUN=quantile, X=qb.lb.ub.optimal.mat[, 1:doses.length], MARGIN=2, probs=alpha.low, na.rm=TRUE)
  qb.CI.ub.vec <- apply(FUN=quantile, X=qb.lb.ub.optimal.mat[, (doses.length+1):(doses.length*2)], MARGIN=2, probs=alpha.high, na.rm=TRUE)
  
  end.time.qb1 <- Sys.time()
  
  start.time.jesson1 <- Sys.time()
  
  # Store Jesson et al. PEIs
  jesson.PEI.lb.vec <- qb.jesson.pei$jesson.PEI.per.dose.lb
  jesson.PEI.ub.vec <- qb.jesson.pei$jesson.PEI.per.dose.ub
  
  # Compute the CIs
  # For Jesson et al.
  jesson.CI.lb.vec <- apply(FUN=quantile, X=jesson.boot.lb.mat, MARGIN=2, probs=alpha.low)
  jesson.CI.ub.vec <- apply(FUN=quantile, X=jesson.boot.ub.mat, MARGIN=2, probs=alpha.high)
  
  end.time.jesson1 <- Sys.time()
  
  # Stop measuring execution time
  exec.time <- toc()
  
  tot.exec.time.qb <- exec.time.qb.PEI + sum(exec.times.qb.vec) + as.numeric(end.time.qb1 - start.time.qb1, units="secs")
  tot.exec.time.jesson <- exec.time.jesson.PEI + sum(exec.times.jesson.vec) + as.numeric(end.time.jesson1 - start.time.jesson1, units="secs")
  
  
  if (data.name == "cmr") {
    
    # Rescale results
    doses.rescaled <- doses * attr(scaled.t, 'scaled:scale') + attr(scaled.t, 'scaled:center')
    
    qb.PEI.lb.vec.rescaled <- qb.PEI.lb.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    qb.PEI.ub.vec.rescaled <- qb.PEI.ub.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    qb.unconf.vec.rescaled <- qb.unconf.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    jesson.PEI.lb.vec.rescaled <- jesson.PEI.lb.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    jesson.PEI.ub.vec.rescaled <- jesson.PEI.ub.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    jesson.unconf.vec.rescaled <- jesson.unconf.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    
    qb.CI.lb.vec.rescaled <- qb.CI.lb.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    qb.CI.ub.vec.rescaled <- qb.CI.ub.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    jesson.CI.lb.vec.rescaled <- jesson.CI.lb.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    jesson.CI.ub.vec.rescaled <- jesson.CI.ub.vec * attr(scaled.Y, 'scaled:scale') + attr(scaled.Y, 'scaled:center')
    
    # Average lengths
    
    # Continuous QB
    message("Continuous QB")
    qb.PEI.avg.length <- mean(qb.PEI.ub.vec.rescaled - qb.PEI.lb.vec.rescaled)
    qb.CI.avg.length <- mean(qb.CI.ub.vec.rescaled - qb.CI.lb.vec.rescaled)
    
    print(paste("qb.PEI.avg.length =", round(qb.PEI.avg.length, 3)))
    print(paste("qb.CI.avg.length =", round(qb.CI.avg.length, 3)))
    
    # Jesson et al.
    message("Jesson et al.")
    jesson.PEI.avg.length <- mean(jesson.PEI.ub.vec.rescaled - jesson.PEI.lb.vec.rescaled)
    jesson.CI.avg.length <- mean(jesson.CI.ub.vec.rescaled - jesson.CI.lb.vec.rescaled)
    
    print(paste("jesson.PEI.avg.length =", round(jesson.PEI.avg.length, 3)))
    print(paste("jesson.CI.avg.length =", round(jesson.CI.avg.length, 3)))
    
    if (use.stabilization) {
      cont.qb.method <- "Continuous SAIPW-QB"
    } else {
      cont.qb.method <- "Continuous AIPW-QB"
    }
    
    qb.bounds.df <- rbind(data.frame(lb=qb.PEI.lb.vec.rescaled,
                                     ub=qb.PEI.ub.vec.rescaled,
                                     unconf=qb.unconf.vec.rescaled,
                                     apo=NA,
                                     dose=doses.rescaled,
                                     name=rep("Cont.QB.PEI", doses.length),
                                     method=rep(cont.qb.method, doses.length),
                                     type=rep("PEI", doses.length)),
                          data.frame(lb=qb.CI.lb.vec.rescaled,
                                     ub=qb.CI.ub.vec.rescaled,
                                     unconf=qb.unconf.vec.rescaled,
                                     apo=NA,
                                     dose=doses.rescaled,
                                     name=rep("Cont.QB.CI", doses.length),
                                     method=rep(cont.qb.method, doses.length),
                                     type=rep("CI", doses.length)))
    
    jesson.bounds.df <- rbind(data.frame(lb=jesson.PEI.lb.vec.rescaled,
                                         ub=jesson.PEI.ub.vec.rescaled,
                                         unconf=jesson.unconf.vec.rescaled,
                                         apo=NA,
                                         dose=doses.rescaled,
                                         name=rep("Jesson.PEI", doses.length),
                                         method=rep("Jesson et al.", doses.length),
                                         type=rep("PEI", doses.length)),
                              data.frame(lb=jesson.CI.lb.vec.rescaled,
                                         ub=jesson.CI.ub.vec.rescaled,
                                         unconf=jesson.unconf.vec.rescaled,
                                         apo=NA,
                                         dose=doses.rescaled,
                                         name=rep("Jesson.CI", doses.length),
                                         method=rep("Jesson et al.", doses.length),
                                         type=rep("CI", doses.length)))
    
  } else if (data.name == "simul") {
    
    apo <- rep(0, length(doses))
    
    for (k in 1:length(doses)) {
      apo[k] <- apo_t(doses[k], zeta, gamma_X, gamma_U, cov_XU, Sigma_X)
    }
    
    # Average lengths
    
    # Continuous QB
    message("Continuous QB")
    qb.PEI.avg.length <- mean(qb.PEI.ub.vec - qb.PEI.lb.vec)
    qb.CI.avg.length <- mean(qb.CI.ub.vec - qb.CI.lb.vec)
    
    print(paste("qb.PEI.avg.length =", round(qb.PEI.avg.length, 3)))
    print(paste("qb.CI.avg.length =", round(qb.CI.avg.length, 3)))
    
    # Jesson et al.
    message("Jesson et al.")
    jesson.PEI.avg.length <- mean(jesson.PEI.ub.vec - jesson.PEI.lb.vec)
    jesson.CI.avg.length <- mean(jesson.CI.ub.vec - jesson.CI.lb.vec)
    
    print(paste("jesson.PEI.avg.length =", round(jesson.PEI.avg.length, 3)))
    print(paste("jesson.CI.avg.length =", round(jesson.CI.avg.length, 3)))
    
    
    apo <- rep(0, length(doses))
    
    for (k in 1:length(doses)) {
      apo[k] <- apo_t(doses[k], zeta, gamma_X, gamma_U, cov_XU, Sigma_X)
    }
    
    if (use.stabilization) {
      cont.qb.method <- "Continuous SAIPW-QB"
    } else {
      cont.qb.method <- "Continuous AIPW-QB"
    }
    
    qb.bounds.df <- rbind(data.frame(lb=qb.PEI.lb.vec,
                                     ub=qb.PEI.ub.vec,
                                     unconf=qb.unconf.vec,
                                     apo=apo_for_plot_fun(doses),
                                     dose=doses,
                                     name=rep("Cont.QB.PEI", doses.length),
                                     method=rep(cont.qb.method, doses.length),
                                     type=rep("PEI", doses.length)),
                          data.frame(lb=qb.CI.lb.vec,
                                     ub=qb.CI.ub.vec,
                                     unconf=qb.unconf.vec,
                                     apo=apo_for_plot_fun(doses),
                                     dose=doses,
                                     name=rep("Cont.QB.CI", doses.length),
                                     method=rep(cont.qb.method, doses.length),
                                     type=rep("CI", doses.length)))
    
    jesson.bounds.df <- rbind(data.frame(lb=jesson.PEI.lb.vec,
                                         ub=jesson.PEI.ub.vec,
                                         unconf=jesson.unconf.vec,
                                         apo=apo_for_plot_fun(doses),
                                         dose=doses,
                                         name=rep("Jesson.PEI", doses.length),
                                         method=rep("Jesson et al.", doses.length),
                                         type=rep("PEI", doses.length)),
                              data.frame(lb=jesson.CI.lb.vec,
                                         ub=jesson.CI.ub.vec,
                                         unconf=jesson.unconf.vec,
                                         apo=apo_for_plot_fun(doses),
                                         dose=doses,
                                         name=rep("Jesson.CI", doses.length),
                                         method=rep("Jesson et al.", doses.length),
                                         type=rep("CI", doses.length)))
    
  } else {
    stop("data.name must be 'cmr' or 'simul'")
  }
  
  
  # Store the results on the Monte-Carlo sample
  qb.mc.bounds[[mc.ind]] <- qb.bounds.df
  jesson.mc.bounds[[mc.ind]] <- jesson.bounds.df
  
  qb.exec.times[mc.ind] <- tot.exec.time.qb
  jesson.exec.times[mc.ind] <- tot.exec.time.jesson
  
  # Save data in a file
  version <- "temp"
  
  # Save the bounds
  qb.file.name <- paste("./results/qb", data.name, "APO_bounds", "B", B, version, sep="_")
  saveRDS(qb.mc.bounds, file=paste0(qb.file.name, ".RData"))
  # Save the execution times
  qb.exec.time.file.name <- paste("./results/qb", data.name, "exec_time", "B", B, version, sep="_")
  saveRDS(qb.exec.times, file=paste0(qb.exec.time.file.name, ".RData"))
  
  # Save the bounds
  jesson.file.name <- paste("./results/jesson", data.name, "APO_bounds", "B", B, version, sep="_")
  saveRDS(jesson.mc.bounds, file=paste0(jesson.file.name, ".RData"))
  # Save the execution times
  jesson.exec.time.file.name <- paste("./results/jesson", data.name, "exec_time", "B", B, version, sep="_")
  saveRDS(jesson.exec.times, file=paste0(jesson.exec.time.file.name, ".RData"))
  
}


if (parallel.computation.jesson) {
  parallel::stopCluster(cl)  # Stop parallel computation
}


# Compute the mean coverage for each Monte-Carlo sample
qb.coverage.vec <- rep(NA, n.MC)
jesson.coverage.vec <- rep(NA, n.MC)
# the mean CI length
qb.CI.length.vec <- rep(NA, n.MC)
jesson.CI.length.vec <- rep(NA, n.MC)
# and the mean PEI length
qb.PEI.length.vec <- rep(NA, n.MC)
jesson.PEI.length.vec <- rep(NA, n.MC)

other.qb.df <- NULL
other.jesson.df <- NULL

for (i in 1:n.MC) {
  
  message(paste("Monte-Carlo sample:", i, "/", n.MC))
  
  # For continuous QB
  qb.mc.bounds.i <- qb.mc.bounds[[i]]
  qb.mc.bounds.i.ci <- qb.mc.bounds.i[qb.mc.bounds.i$type == "CI", ]
  qb.mc.bounds.i.pei <- qb.mc.bounds.i[qb.mc.bounds.i$type == "PEI", ]
  
  qb.coverage.vec[i] <- mean(qb.mc.bounds.i.ci$lb <= qb.mc.bounds.i.ci$apo & qb.mc.bounds.i.ci$apo <= qb.mc.bounds.i.ci$ub)
  qb.CI.length.vec[i] <- mean(qb.mc.bounds.i.ci$ub - qb.mc.bounds.i.ci$lb)
  qb.PEI.length.vec[i] <- mean(qb.mc.bounds.i.pei$ub - qb.mc.bounds.i.pei$lb)
  
  qb.mc.bounds.i$MC.ind <- i
  other.qb.df <- rbind(other.qb.df, qb.mc.bounds.i)
  
  # For Jesson
  jesson.mc.bounds.i <- jesson.mc.bounds[[i]]
  jesson.mc.bounds.i.ci <- jesson.mc.bounds.i[jesson.mc.bounds.i$type == "CI", ]
  jesson.mc.bounds.i.pei <- jesson.mc.bounds.i[jesson.mc.bounds.i$type == "PEI", ]
  
  jesson.coverage.vec[i] <- mean(jesson.mc.bounds.i.ci$lb <= jesson.mc.bounds.i.ci$apo & jesson.mc.bounds.i.ci$apo <= jesson.mc.bounds.i.ci$ub)
  jesson.CI.length.vec[i] <- mean(jesson.mc.bounds.i.ci$ub - jesson.mc.bounds.i.ci$lb)
  jesson.PEI.length.vec[i] <- mean(jesson.mc.bounds.i.pei$ub - jesson.mc.bounds.i.pei$lb)
  
  jesson.mc.bounds.i$MC.ind <- i
  other.jesson.df <- rbind(other.jesson.df, jesson.mc.bounds.i)
}


# Wilcoxon signed-rank test
wilcox.test.ci.length <- wilcox.test(x=qb.CI.length.vec, y=jesson.CI.length.vec, alternative="two.sided", mu=0, paired=TRUE)
wilcox.test.ci.length

# Wilcoxon signed-rank test
wilcox.test.pei.length <- wilcox.test(x=qb.PEI.length.vec, y=jesson.PEI.length.vec, alternative="two.sided", mu=0, paired=TRUE)
wilcox.test.pei.length

qb.PEI.sd <- sd(qb.PEI.length.vec)
jesson.PEI.sd <- sd(jesson.PEI.length.vec)

paste("QB PEI standard deviation =", qb.PEI.sd)
paste("Jesson PEI standard deviation =", jesson.PEI.sd)


# Save data in a file
version <- paste0("v", version.arg)

# Save the bounds
qb.file.name <- paste("./results/qb", data.name, "APO_bounds", "B", B, version, sep="_")
saveRDS(qb.mc.bounds, file=paste0(qb.file.name, ".RData"))
# Save the execution times
qb.exec.time.file.name <- paste("./results/qb", data.name, "exec_time", "B", B, version, sep="_")
saveRDS(qb.exec.times, file=paste0(qb.exec.time.file.name, ".RData"))
# Save the coverage
qb.cov.file.name <- paste("./results/qb", data.name, "cov", "B", B, version, sep="_")
saveRDS(qb.coverage.vec, file=paste0(qb.cov.file.name, ".RData"))

# Save the bounds
jesson.file.name <- paste("./results/jesson", data.name, "APO_bounds", "B", B, version, sep="_")
saveRDS(jesson.mc.bounds, file=paste0(jesson.file.name, ".RData"))
# Save the execution times
jesson.exec.time.file.name <- paste("./results/jesson", data.name, "exec_time", "B", B, version, sep="_")
saveRDS(jesson.exec.times, file=paste0(jesson.exec.time.file.name, ".RData"))
# Save the coverage
jesson.cov.file.name <- paste("./results/jesson", data.name, "cov", "B", B, version, sep="_")
saveRDS(jesson.coverage.vec, file=paste0(jesson.cov.file.name, ".RData"))

# Save used gammas
gammas.file.name <- paste("./results/gammas", data.name, version, sep="_")
saveRDS(gammas, file=paste0(gammas.file.name, ".RData"))
