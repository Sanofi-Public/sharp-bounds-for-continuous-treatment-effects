# Functions used particularly by Jesson et al.'s method


# Heaviside function
heaviside.fun <- function(x) {
  return(ifelse(x >= 0, 1, 0))
}


# Function to compute the PEI
APO.Jesson.PEI <- function(dose, lambda, capo.method, nn.init=NULL, d1.d2.ind, unconf.capo.model=NULL,
                           X.tensor, X.sample, Y.sample.len, use.parallel=TRUE, device, verbose=TRUE) {
  
  n.obs <- dim(X.tensor)[1]
  
  message("Estimating the PEI")
  message("Sampling outcomes according to p(y|x,t)")
  
  Y.samples <- matrix(NA, nrow=n.obs, ncol=Y.sample.len)
  
  # The estimated density for particular values of X and t
  Y.density1 <- nn.init$resp.gmm2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                  treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data1.ind)), device=device), dim=1))
  Y.samples1 <- Y.density1$sample(Y.sample.len)
  Y.samples1 <- t(as.matrix(Y.samples1))  # The samples are now on each row (n.data1 * Y.sample.len)

  # The estimated density for particular values of X and t
  Y.density2 <- nn.init$resp.gmm1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                  treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data2.ind)), device=device), dim=1))
  Y.samples2 <- Y.density2$sample(Y.sample.len)
  Y.samples2 <- t(as.matrix(Y.samples2))  # The samples are now on each row (n.data2 * Y.sample.len)

  # Reorganize the results
  Y.samples[d1.d2.ind$data1.ind, ] <- Y.samples1
  Y.samples[d1.d2.ind$data2.ind, ] <- Y.samples2
  
  
  message("Estimating the CAPOs under ignorability")
  
  # If the method is "Neural Network Gaussian Mixture Model"
  if (capo.method == "nn_gmm") {
    
    mu.tilde <- rep(NA, n.obs)
    
    mix.prob1 <- Y.density1$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean1 <- Y.density1$.component_distribution$mean  # The fitted means for each component
    mu.tilde1 <- torch_sum(mix.prob1 * dist.mean1, dim=-1)  # Compute the CAPO under ignorability assumption
    mu.tilde1 <- as.array(mu.tilde1)

    mix.prob2 <- Y.density2$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean2 <- Y.density2$.component_distribution$mean  # The fitted means for each component
    mu.tilde2 <- torch_sum(mix.prob2 * dist.mean2, dim=-1)  # Compute the CAPO under ignorability assumption
    mu.tilde2 <- as.array(mu.tilde2)

    # Reorganize the predictions
    mu.tilde[d1.d2.ind$data1.ind] <- c(mu.tilde1)
    mu.tilde[d1.d2.ind$data2.ind] <- c(mu.tilde2)
    
    # If the method is "Regression forest"
  } else if (capo.method == "regression_forest") {
    
    # Compute mu.tilde, the CAPO with unconfoundedness assumption
    mu.tilde <- predict(unconf.capo.model, data.frame(X=X.sample, t=rep(dose, X.sample.len)))$predictions
    
  } else {
    stop('capo.method must be "nn_gmm" or "regression_forest"')
  }
  
  unconf.APO.estimate <- mean(mu.tilde)  # Take the mean of the estimated CAPOs to get the APO estimate under ignorability assumption for a given dose
  
  
  message("Estimating the CAPO bounds")

  if (use.parallel) {  # If we want to use parallel computation
    
    if (verbose) {
      # Initialize another progress bar
      pb <- txtProgressBar(max=n.obs, style=3)
      progress <- function(n) {
        setTxtProgressBar(pb, n)
        if (n == n.obs) {close(pb)}}
      opts <- list(progress=progress)
    } else {
      opts <- NULL
    }
    
    # Find y.down and y.up thanks to grid search (Algorithm 1)
    CAPO.estimates <- foreach(k=1:n.obs, .combine="rbind",
                              .options.snow=opts) %dopar% {
                                
                                CAPO.bounds <- grid.search.optimizer(lambda=lambda,
                                                                     Y.sample=Y.samples[k, ],
                                                                     mu.tilde=mu.tilde[k])
                                
                                CAPO.bounds
                              }
    
  } else {
    
    if (verbose) {
      pb <- txtProgressBar(max=n.obs, style=3)
    }
    
    # Find y.down and y.up thanks to grid search (Algorithm 1)
    CAPO.estimates <- foreach(k=1:n.obs, .combine="rbind") %do% {
                                
                                CAPO.bounds <- grid.search.optimizer(lambda=lambda,
                                                                     Y.sample=Y.samples[k, ],
                                                                     mu.tilde=mu.tilde[k])
                                
                                if (verbose) {
                                  setTxtProgressBar(pb, k)
                                }
                                
                                CAPO.bounds
                              }
    
    if (verbose) {
      close(pb)
    }
  }
  
  message("Computing the APO bounds")
  # Compute the lower and upper bounds of the PEI for the APO
  APO.lower <- mean(CAPO.estimates[, 3])
  APO.upper <- mean(CAPO.estimates[, 4])
  
  return(list(PEI.bounds=c(APO.lower, APO.upper), unconf.APO.estimate=unconf.APO.estimate))
}


# Function to estimate the CAPO bounds
CAPO.estim <- function(y.vec, y.H, optim, mu.tilde, lambda, w.fun) {
  # If we want to solve the minimization problem
  if (optim == "min") {
    
    w.y <- w.fun(rep(y.H, length(y.vec)) - y.vec)
    capo <- mu.tilde + mean(w.y * (y.vec - rep(mu.tilde, length(y.vec)))) / (1/(lambda**2 - 1) + mean(w.y))
    
  } else if (optim == "max") {  # Else, if we want to solve the maximization problem
    
    w.y <- w.fun(y.vec - rep(y.H, length(y.vec)))
    capo <- mu.tilde + mean(w.y * (y.vec - rep(mu.tilde, length(y.vec)))) / (1/(lambda**2 - 1) + mean(w.y))
    
  } else {
    stop("optim must be 'min' or 'max'")
  }
  
  return(capo)
}


# Find values of y.down and y.up thanks to Algorithm 1 from Jesson et al.
# lambda is the sensitivity parameter
# Y.sample is the sample obtained for a given couple (X, dose)
# mu.tilde is the CAPO estimated under ignorability assumption for a given couple (X, dose)
grid.search.optimizer <- function(lambda, Y.sample, mu.tilde) {
  
  CAPO.up.fun <- function(y) CAPO.estim(y.vec=Y.sample, y.H=y, optim="max", mu.tilde=mu.tilde, lambda=lambda, w.fun=heaviside.fun)
  CAPO.down.fun <- function(y) CAPO.estim(y.vec=Y.sample, y.H=y, optim="min", mu.tilde=mu.tilde, lambda=lambda, w.fun=heaviside.fun)

  kappa.up.vec <- sapply(X=Y.sample, FUN=CAPO.up.fun)
  kappa.down.vec <- sapply(X=Y.sample, FUN=CAPO.down.fun)

  return(c(0, 0, min(kappa.down.vec), max(kappa.up.vec)))
  
  mu.up <- -Inf
  mu.down <- Inf
  
  y.up <- 0
  y.down <- 0
  
  Y.sample <- unique(Y.sample)  # Make sure we only have unique values
  Y.sample.len <- length(Y.sample)
  
  for (i in 1:Y.sample.len) {
    Y.i <- Y.sample[i]
    
    kappa.up <- CAPO.estim(y.vec=Y.sample, y.H=Y.i, optim="max", mu.tilde=mu.tilde, lambda=lambda, w.fun=heaviside.fun)
    kappa.down <- CAPO.estim(y.vec=Y.sample, y.H=Y.i, optim="min", mu.tilde=mu.tilde, lambda=lambda, w.fun=heaviside.fun)
    
    if (kappa.up > mu.up) {
      mu.up <- kappa.up
      y.up <- Y.i
    }
    
    if (kappa.down < mu.down) {
      mu.down <- kappa.down
      y.down <- Y.i
    }
  }
  
  return(c(y.down, y.up, mu.down, mu.up))
}