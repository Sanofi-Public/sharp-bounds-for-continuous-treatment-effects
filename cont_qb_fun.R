# Functions used particularly by our method


# Epanechnikov kernel
epa.kernel <- function(u) {
  return(0.75 * (1 - u**2) * (abs(u) <= 1))
}


# Kernel with window h
K.h <- function(u, h) {
  return(epa.kernel(u/h) / h)
}


# Function to compute conditional quantiles
cond.quant.estim <- function(X, t, resp.gmm, tau, device, verbose) {
  
  X.tensor <- X
  t.tensor <- t
  
  message("Estimating the conditional quantiles using p(y|x,t)")

  if (verbose) {
    pb <- txtProgressBar(max=dim(X.tensor)[1], style=3)
  }
  
  Q.predict <- foreach(i=1:dim(X.tensor)[1], .combine="rbind") %do% {
    
    # I want to plot the cdf of y|x,t for individual ind in data2
    resp.dens <- resp.gmm(covariates=X.tensor[i, ], treatment=t.tensor[i])
    
    min.x <- -10
    max.x <- 10
    
    cdf.min.one_tau.fun <- function(x) as.numeric(resp.dens$cdf(torch_tensor(x, device=device))$to(device="cpu"))-(1-tau)
    cdf.min.tau.fun <- function(x) as.numeric(resp.dens$cdf(torch_tensor(x, device=device))$to(device="cpu"))-tau
    
    one_tau.quant <- uniroot(f=cdf.min.one_tau.fun, interval=c(min.x, max.x))$root
    tau.quant <- uniroot(f=cdf.min.tau.fun, interval=c(min.x, max.x))$root
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    
    c(one_tau.quant, tau.quant)
  }
  
  if (verbose) {
    close(pb)
  }
  
  return(Q.predict)
}


# Function to compute the PEI for each provided window h
APO.QB.PEI.per.window <- function(X, Y, t, dose, alpha=0.05, lambda=1,
                                  window=0.1, windows.range=NULL, B.param=50,
                                  cond.quant.method=c("quantile_forest"),
                                  boot.training=FALSE, max.iter=1500,
                                  optimal.params.gps.nn=list(K.optim=25, lr.optim=0.0015, dim.hidden.optim=32),
                                  optimal.params.capo.nn=list(K.optim=25, lr.optim=0.0015, dim.hidden.optim=32),
                                  stabilization=TRUE,
                                  use.AIPW=FALSE, AIPW.outcome.algo=c("nn_GMM"),
                                  cond.dens.method=c("nn_GMM"), nn.init=NULL, d1.d2.ind=NULL, h1=0.2, h2=0.2, K=NULL, K.range=1:20, K.nrep=3,
                                  xi.method=c("regression_forest"), xi.models,
                                  compute.CI=TRUE, B=1000, use.parallel=TRUE, e.variables=NULL,
                                  use.cross.fitting=TRUE, nb.folds=5, epsilon=10**(-12), seed=1, device, verbose=TRUE) {
  
  # Put data in a dataframe
  simul.data <- data.frame(X, Y, t)
  
  # Number of observations
  n.obs <- nrow(simul.data)
  
  tau <- lambda / (1 + lambda)
  windows.len <- length(windows.range)  # Number of windows that we want to test
  
  # If windows.range was provided
  if (!is.null(windows.range)) {
    
    PEI.bounds.per.window <- compute.PEI.cont.QB(dataset=simul.data,
                                                 dose=dose,
                                                 lambda=lambda,
                                                 windows=windows.range,
                                                 stabilization=stabilization,
                                                 use.AIPW=use.AIPW,
                                                 AIPW.outcome.algo=AIPW.outcome.algo,
                                                 cond.dens.method=cond.dens.method,
                                                 xi.method=xi.method,
                                                 xi.models=xi.models,
                                                 nn.init=nn.init,
                                                 d1.d2.ind=d1.d2.ind,
                                                 h1=h1, h2=h2, K=K,
                                                 e.variables=e.variables,
                                                 use.cross.fitting=use.cross.fitting,
                                                 nb.folds=nb.folds,
                                                 epsilon=epsilon, device=device)
    
    PEI.min.max <- NULL
    
  }
  
  return(list(PEI.bounds.per.window=PEI.bounds.per.window, PEI.min.max=PEI.min.max,
              dose=dose, use.AIPW=use.AIPW, AIPW.outcome.algo=AIPW.outcome.algo,
              cond.dens.method=cond.dens.method, h1=h1, h2=h2, K=K,
              use.cross.fitting=use.cross.fitting, nb.folds=nb.folds, e.variables=e.variables,
              windows.range=windows.range, windows.len=windows.len,
              simul.data=simul.data, n.obs=n.obs, lambda=lambda, tau=tau, epsilon=epsilon))
}


compute.PEI.cont.QB <- function(dataset, dose, lambda, windows,
                                stabilization,
                                use.AIPW, AIPW.outcome.algo,
                                cond.dens.method, nn.init,
                                xi.method, xi.models,
                                d1.d2.ind,
                                h1=0.2, h2=0.2, K=3,
                                e.variables,
                                use.cross.fitting, nb.folds=5,
                                epsilon=10**(-6), device) {
  
  n.data <- nrow(dataset)
  windows.length <- length(windows)
  
  X <- subset(dataset, select=-c(Y, t))
  n.cov <- ncol(X)
  
  # Reorganize the predictions for the conditional quantiles
  Q.predict <- matrix(NA, nrow=n.data, ncol=2)
  Q.predict[d1.d2.ind$data1.ind, ] <- nn.init$Q.predict1
  Q.predict[d1.d2.ind$data2.ind, ] <- nn.init$Q.predict2
  
  
  ################ Generalized Propensity Score estimation (p(t|x)) ################
  
  if (cond.dens.method == "true") {
    
    # Use the true GPS (only for simulated data from Jesson et al.)
    cond.density <- true.cond.dens(dataset$t*100, dataset$X)  # Use t*100 to transform into integers between 1 and 100
    
  } else if (cond.dens.method == "nn_GMM") {
    
    # GPS estimation via Gaussian Mixture Model fitted via neural networks
    
    X.tensor <- torch_tensor(as.matrix(X), device=device)
    t.tensor <- torch_tensor(dataset$t, device=device)
    
    cond.density <- rep(NA, n.data)

    # On D1
    cond.density1 <- as.array(exp(nn.init$gps.gmm2(covariates=X.tensor[d1.d2.ind$data1.ind, ])$log_prob(x=t.tensor[d1.d2.ind$data1.ind])$to(device="cpu")))
    # On D2
    cond.density2 <- as.array(exp(nn.init$gps.gmm1(covariates=X.tensor[d1.d2.ind$data2.ind, ])$log_prob(x=t.tensor[d1.d2.ind$data2.ind])$to(device="cpu")))
    
    # Reorganize the predictions
    cond.density[d1.d2.ind$data1.ind] <- cond.density1
    cond.density[d1.d2.ind$data2.ind] <- cond.density2
    
  } else if (cond.dens.method == "kernel") {
    
    # Use kernel density estimation (only for 1-dimensional covariates)
    
    if (use.cross.fitting) {
      
      # Cross-fitting
      nb.data.per.fold <- n.data %/% nb.folds
      cond.density <- rep(NA, n.data)
      
      for (j in 1:nb.folds) {
        
        if (j == nb.folds) {  # For the last fold, we want to be sure we have used all the data
          prediction.fold.ind <- (1+nb.data.per.fold*(j-1)):n.data
          prediction.fold <- dataset[prediction.fold.ind, ]
          training.fold <- dataset[-prediction.fold.ind, ]
        } else {
          prediction.fold.ind <- (1+nb.data.per.fold*(j-1)):(nb.data.per.fold*j)
          prediction.fold <- dataset[prediction.fold.ind, ]
          training.fold <- dataset[-prediction.fold.ind, ]
        }
        
        # Nadaraya-Watson-like conditional density estimator
        num.fun <- function(x) sum(K.h(training.fold$t - dose, h1) * K.h(training.fold$X - x, h2))
        denom.fun <- function(x) sum(K.h(training.fold$X - x, h2))
        
        num <- vapply(X=prediction.fold$X, FUN=num.fun, FUN.VALUE=as.numeric(1))
        denom <- vapply(X=prediction.fold$X, FUN=denom.fun, FUN.VALUE=as.numeric(1))
        
        cond.density[prediction.fold.ind] <- num / denom + epsilon  # epsilon = 10**(-12) added not to have 0
        
      }
      
    } else {
      
      # Nadaraya-Watson-like conditional density estimator
      num.fun <- function(x) sum(K.h(dataset$t - dose, h1) * K.h(dataset$X - x, h2))
      denom.fun <- function(x) sum(K.h(dataset$X - x, h2))
      
      num <- vapply(X=dataset$X, FUN=num.fun, FUN.VALUE=as.numeric(1))
      denom <- vapply(X=dataset$X, FUN=denom.fun, FUN.VALUE=as.numeric(1))
      
      cond.density <- num / denom + epsilon
      
    }
    
  } else if (cond.dens.method == "flexmix_GMM") {
    
    # GPS estimation via mixture of linear regressions
    
    if (use.cross.fitting) {
      
      # Cross-fitting
      nb.data.per.fold <- n.data %/% nb.folds
      cond.density <- rep(NA, n.data)  # Initialize conditional density vector
      
      for (j in 1:nb.folds) {
        
        if (j == nb.folds) {  # For the last fold, we want to be sure we have used all the data
          prediction.fold.ind <- (1+nb.data.per.fold*(j-1)):n.data
          prediction.fold <- dataset[prediction.fold.ind, ]
          training.fold <- dataset[-prediction.fold.ind, ]
        } else {
          prediction.fold.ind <- (1+nb.data.per.fold*(j-1)):(nb.data.per.fold*j)
          prediction.fold <- dataset[prediction.fold.ind, ]
          training.fold <- dataset[-prediction.fold.ind, ]
        }
        
        # Estimate the parameters of the Gaussian Mixture Model (GMM)
        flexmix.gmm <- flexmix(formula=t ~ . - Y, data=training.fold, k=K)
        
        # Predict the values of the generalized propensity score
        cond.density[prediction.fold.ind] <- est.dens.mat(Y=prediction.fold$t,
                                                          X=subset(prediction.fold, select=-c(Y, t)),
                                                          flexmix.model=flexmix.gmm)
        
      }
      
    } else {
      
      # Create the formula that will be given to flexmix()
      # Estimate the parameters of the Gaussian Mixture Model (GMM)
      flexmix.gmm <- flexmix(formula=t ~ . - Y, data=dataset, k=K)
      # Predict the values of the generalized propensity score
      # cond.density <- est.dens.mat(Y=dose, X=subset(dataset.ordered, select=-c(Y, t)), flexmix.model=flexmix.gmm)
      cond.density <- est.dens.mat(Y=dataset$t, X=X, flexmix.model=flexmix.gmm)
    }
    
  } else {
    
    stop("cond.dens.method must be 'true', 'kernel', 'flexmix_GMM' or 'nn_GMM'")
    
  }
  
  # Weight trimming
  trim.lq <- 0.1
  trim.quant <- quantile(cond.density, probs=c(trim.lq))
  # Replace all weights with the quantiles of order trim.lq if they are lower than trim.lq
  cond.density[cond.density < trim.quant[1]] <- trim.quant[1]

  ################ CAPO estimation (E[Y|X=x,T=t]) under ignorability assumption ################
  
  if (use.AIPW) {
    if (AIPW.outcome.algo == "lm") {
      
      # Fit the outcome model
      outcome.model <- lm(Y ~ ., dataset)
      
      # Predict the outcome for the dose assigned to each individual
      Y.predict.all <- predict(outcome.model, dataset)
      
      # Predict the outcome for a specific dose
      predict.dataset <- dataset
      predict.dataset$t <- dose
      Y.predict.dose <- predict(outcome.model, predict.dataset)
      
    } else if (AIPW.outcome.algo == "regression_forest") {
      
      # Fit the outcome model
      outcome.model <- regression_forest(X=data.frame(X=X, t=dataset$t), Y=dataset$Y, num.trees=100)  # To be changed if X is not univariate
      
      # Predict the outcome for the dose assigned to each individual
      Y.predict.all <- predict(outcome.model, data.frame(X, t=dataset$t))$predictions
      
      # Predict the outcome for a specific dose
      Y.predict.dose <- predict(outcome.model, data.frame(X, t=rep(dose, nrow(dataset))))$predictions
      
    } else if (AIPW.outcome.algo == "nn_GMM") {
      
      X.tensor <- torch_tensor(as.matrix(X), device=device)
      t.tensor <- torch_tensor(dataset$t, device=device)
      
      Y.predict.all <- rep(NA, n.data)
      Y.predict.dose <- rep(NA, n.data)
      
      
      # The estimated density for particular values of X and t
      Y.density.all1 <- nn.init$resp.gmm2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                          treatment=t.tensor[d1.d2.ind$data1.ind])
      Y.density.all2 <- nn.init$resp.gmm1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                          treatment=t.tensor[d1.d2.ind$data2.ind])
      
      mix.prob.all1 <- Y.density.all1$.mixture_distribution$probs  # The fitted probabilities of each component
      dist.mean.all1 <- Y.density.all1$.component_distribution$mean  # The fitted means for each component
      Y.predict.all1 <- torch_sum(mix.prob.all1 * dist.mean.all1, dim=-1)
      Y.predict.all1 <- as.array(Y.predict.all1$to(device="cpu"))
      
      mix.prob.all2 <- Y.density.all2$.mixture_distribution$probs  # The fitted probabilities of each component
      dist.mean.all2 <- Y.density.all2$.component_distribution$mean  # The fitted means for each component
      Y.predict.all2 <- torch_sum(mix.prob.all2 * dist.mean.all2, dim=-1)
      Y.predict.all2 <- as.array(Y.predict.all2$to(device="cpu"))
      
      # Reorganize the predictions
      Y.predict.all[d1.d2.ind$data1.ind] <- c(Y.predict.all1)
      Y.predict.all[d1.d2.ind$data2.ind] <- c(Y.predict.all2)
      
      
      # The estimated density for particular values of X and t
      Y.density.dose1 <- nn.init$resp.gmm2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                           treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data1.ind)), device=device), dim=1))
      Y.density.dose2 <- nn.init$resp.gmm1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                           treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data2.ind)), device=device), dim=1))
      
      mix.prob.dose1 <- Y.density.dose1$.mixture_distribution$probs  # The fitted probabilities of each component
      dist.mean.dose1 <- Y.density.dose1$.component_distribution$mean  # The fitted means for each component
      Y.predict.dose1 <- torch_sum(mix.prob.dose1 * dist.mean.dose1, dim=-1)
      Y.predict.dose1 <- as.array(Y.predict.dose1$to(device="cpu"))
      
      mix.prob.dose2 <- Y.density.dose2$.mixture_distribution$probs  # The fitted probabilities of each component
      dist.mean.dose2 <- Y.density.dose2$.component_distribution$mean  # The fitted means for each component
      Y.predict.dose2 <- torch_sum(mix.prob.dose2 * dist.mean.dose2, dim=-1)
      Y.predict.dose2 <- as.array(Y.predict.dose2$to(device="cpu"))
      
      # Reorganize the predictions
      Y.predict.dose[d1.d2.ind$data1.ind] <- c(Y.predict.dose1)
      Y.predict.dose[d1.d2.ind$data2.ind] <- c(Y.predict.dose2)
      
    } else {
      stop("AIPW.outcome.algo must be 'lm' or 'regression_forest'")
    }
    
  } else {
    Y.predict.dose <- rep(0, n.data)
    Y.predict.all <- rep(0, n.data)
  }
  
  
  if (is.null(xi.method)) {
    
  } else if (xi.method == "regression_forest") {
    
    # Predictions of the xi function that brings double robustness
    # Model fitting on D1
    X.t.d1 <- as.matrix(subset(dataset, select=-c(Y))[d1.d2.ind$data1.ind, ])
    # Model fitting on D2
    X.t.d2 <- as.matrix(subset(dataset, select=-c(Y))[d1.d2.ind$data2.ind, ])
    
    xi.t.predict.down <- rep(NA, n.data)
    xi.t.predict.up <- rep(NA, n.data)
    # Predict on D1 with model fitted on D2
    xi.t.predict.down[d1.d2.ind$data1.ind] <- predict(xi.models$xi.model.down.d2, X.t.d1)$predictions * attr(xi.models$Y.gamma.down.d2, 'scaled:scale') + attr(xi.models$Y.gamma.down.d2, 'scaled:center')
    xi.t.predict.up[d1.d2.ind$data1.ind] <- predict(xi.models$xi.model.up.d2, X.t.d1)$predictions * attr(xi.models$Y.gamma.up.d2, 'scaled:scale') + attr(xi.models$Y.gamma.up.d2, 'scaled:center')
    # Predict on D2 with model fitted on D1
    xi.t.predict.down[d1.d2.ind$data2.ind] <- predict(xi.models$xi.model.down.d1, X.t.d2)$predictions * attr(xi.models$Y.gamma.down.d1, 'scaled:scale') + attr(xi.models$Y.gamma.down.d1, 'scaled:center')
    xi.t.predict.up[d1.d2.ind$data2.ind] <- predict(xi.models$xi.model.up.d1, X.t.d2)$predictions * attr(xi.models$Y.gamma.up.d1, 'scaled:scale') + attr(xi.models$Y.gamma.up.d1, 'scaled:center')
    
    # Covariates with t replaced with the dose of interest
    X.dose.d1 <- subset(dataset, select=-c(Y))[d1.d2.ind$data1.ind, ]
    X.dose.d1$t <- dose
    X.dose.d1 <- as.matrix(X.dose.d1)
    
    X.dose.d2 <- subset(dataset, select=-c(Y))[d1.d2.ind$data2.ind, ]
    X.dose.d2$t <- dose
    X.dose.d2 <- as.matrix(X.dose.d2)
    
    xi.dose.predict.down <- rep(NA, n.data)
    xi.dose.predict.up <- rep(NA, n.data)
    # Predict on D1 with model fitted on D2
    xi.dose.predict.down[d1.d2.ind$data1.ind] <- predict(xi.models$xi.model.down.d2, X.dose.d1)$predictions * attr(xi.models$Y.gamma.down.d2, 'scaled:scale') + attr(xi.models$Y.gamma.down.d2, 'scaled:center')
    xi.dose.predict.up[d1.d2.ind$data1.ind] <- predict(xi.models$xi.model.up.d2, X.dose.d1)$predictions * attr(xi.models$Y.gamma.up.d2, 'scaled:scale') + attr(xi.models$Y.gamma.up.d2, 'scaled:center')
    # Predict on D2 with model fitted on D1
    xi.dose.predict.down[d1.d2.ind$data2.ind] <- predict(xi.models$xi.model.down.d1, X.dose.d2)$predictions * attr(xi.models$Y.gamma.down.d1, 'scaled:scale') + attr(xi.models$Y.gamma.down.d1, 'scaled:center')
    xi.dose.predict.up[d1.d2.ind$data2.ind] <- predict(xi.models$xi.model.up.d1, X.dose.d2)$predictions * attr(xi.models$Y.gamma.up.d1, 'scaled:scale') + attr(xi.models$Y.gamma.up.d1, 'scaled:center')
    
  } else if (xi.method == "neural_network") {
    
    X.tensor <- torch_tensor(as.matrix(X), device=device)
    t.tensor <- torch_tensor(dataset$t, device=device)
    
    xi.t.predict.down <- rep(NA, n.data)
    xi.t.predict.up <- rep(NA, n.data)
    xi.dose.predict.down <- rep(NA, n.data)
    xi.dose.predict.up <- rep(NA, n.data)
    
    # For the lower bound
    # The estimated density for particular values of X and t
    xi.density.t.down1 <- xi.models$xi.gmm.down2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                                 treatment=t.tensor[d1.d2.ind$data1.ind])
    xi.density.t.down2 <- xi.models$xi.gmm.down1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                                 treatment=t.tensor[d1.d2.ind$data2.ind])
    
    mix.prob.t.down1 <- xi.density.t.down1$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.t.down1 <- xi.density.t.down1$.component_distribution$mean  # The fitted means for each component
    xi.t.predict.down1 <- torch_sum(mix.prob.t.down1 * dist.mean.t.down1, dim=-1)
    xi.t.predict.down1 <- as.array(xi.t.predict.down1$to(device="cpu"))
    
    mix.prob.t.down2 <- xi.density.t.down2$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.t.down2 <- xi.density.t.down2$.component_distribution$mean  # The fitted means for each component
    xi.t.predict.down2 <- torch_sum(mix.prob.t.down2 * dist.mean.t.down2, dim=-1)
    xi.t.predict.down2 <- as.array(xi.t.predict.down2$to(device="cpu"))
    
    # Reorganize the predictions
    xi.t.predict.down[d1.d2.ind$data1.ind] <- c(xi.t.predict.down1) * attr(xi.models$y.train.down2, 'scaled:scale') + attr(xi.models$y.train.down2, 'scaled:center')
    xi.t.predict.down[d1.d2.ind$data2.ind] <- c(xi.t.predict.down2) * attr(xi.models$y.train.down1, 'scaled:scale') + attr(xi.models$y.train.down1, 'scaled:center')
    
    # For the upper bound
    # The estimated density for particular values of X and t
    xi.density.t.up1 <- xi.models$xi.gmm.up2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                             treatment=t.tensor[d1.d2.ind$data1.ind])
    xi.density.t.up2 <- xi.models$xi.gmm.up1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                             treatment=t.tensor[d1.d2.ind$data2.ind])
    
    mix.prob.t.up1 <- xi.density.t.up1$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.t.up1 <- xi.density.t.up1$.component_distribution$mean  # The fitted means for each component
    xi.t.predict.up1 <- torch_sum(mix.prob.t.up1 * dist.mean.t.up1, dim=-1)
    xi.t.predict.up1 <- as.array(xi.t.predict.up1$to(device="cpu"))
    
    mix.prob.t.up2 <- xi.density.t.up2$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.t.up2 <- xi.density.t.up2$.component_distribution$mean  # The fitted means for each component
    xi.t.predict.up2 <- torch_sum(mix.prob.t.up2 * dist.mean.t.up2, dim=-1)
    xi.t.predict.up2 <- as.array(xi.t.predict.up2$to(device="cpu"))
    
    # Reorganize the predictions
    xi.t.predict.up[d1.d2.ind$data1.ind] <- c(xi.t.predict.up1) * attr(xi.models$y.train.up2, 'scaled:scale') + attr(xi.models$y.train.up2, 'scaled:center')
    xi.t.predict.up[d1.d2.ind$data2.ind] <- c(xi.t.predict.up2) * attr(xi.models$y.train.up1, 'scaled:scale') + attr(xi.models$y.train.up1, 'scaled:center')
    
    
    # For a particular value of dose
    # For the lower bound
    # The estimated density for particular values of X and t
    xi.density.dose.down1 <- xi.models$xi.gmm.down2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                                    treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data1.ind)), device=device), dim=1))
    xi.density.dose.down2 <- xi.models$xi.gmm.down1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                                    treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data2.ind)), device=device), dim=1))
    
    mix.prob.dose.down1 <- xi.density.dose.down1$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.dose.down1 <- xi.density.dose.down1$.component_distribution$mean  # The fitted means for each component
    xi.dose.predict.down1 <- torch_sum(mix.prob.dose.down1 * dist.mean.dose.down1, dim=-1)
    xi.dose.predict.down1 <- as.array(xi.dose.predict.down1$to(device="cpu"))
    
    mix.prob.dose.down2 <- xi.density.dose.down2$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.dose.down2 <- xi.density.dose.down2$.component_distribution$mean  # The fitted means for each component
    xi.dose.predict.down2 <- torch_sum(mix.prob.dose.down2 * dist.mean.dose.down2, dim=-1)
    xi.dose.predict.down2 <- as.array(xi.dose.predict.down2$to(device="cpu"))
    
    # Reorganize the predictions
    xi.dose.predict.down[d1.d2.ind$data1.ind] <- c(xi.dose.predict.down1) * attr(xi.models$y.train.down2, 'scaled:scale') + attr(xi.models$y.train.down2, 'scaled:center')
    xi.dose.predict.down[d1.d2.ind$data2.ind] <- c(xi.dose.predict.down2) * attr(xi.models$y.train.down1, 'scaled:scale') + attr(xi.models$y.train.down1, 'scaled:center')
    
    # For the upper bound
    # The estimated density for particular values of X and t
    xi.density.dose.up1 <- xi.models$xi.gmm.up2(covariates=X.tensor[d1.d2.ind$data1.ind, ],
                                                treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data1.ind)), device=device), dim=1))
    xi.density.dose.up2 <- xi.models$xi.gmm.up1(covariates=X.tensor[d1.d2.ind$data2.ind, ],
                                                treatment=torch_squeeze(torch_tensor(rep(dose, length(d1.d2.ind$data2.ind)), device=device), dim=1))
    
    mix.prob.dose.up1 <- xi.density.dose.up1$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.dose.up1 <- xi.density.dose.up1$.component_distribution$mean  # The fitted means for each component
    xi.dose.predict.up1 <- torch_sum(mix.prob.dose.up1 * dist.mean.dose.up1, dim=-1)
    xi.dose.predict.up1 <- as.array(xi.dose.predict.up1$to(device="cpu"))
    
    mix.prob.dose.up2 <- xi.density.dose.up2$.mixture_distribution$probs  # The fitted probabilities of each component
    dist.mean.dose.up2 <- xi.density.dose.up2$.component_distribution$mean  # The fitted means for each component
    xi.dose.predict.up2 <- torch_sum(mix.prob.dose.up2 * dist.mean.dose.up2, dim=-1)
    xi.dose.predict.up2 <- as.array(xi.dose.predict.up2$to(device="cpu"))
    
    # Reorganize the predictions
    xi.dose.predict.up[d1.d2.ind$data1.ind] <- c(xi.dose.predict.up1) * attr(xi.models$y.train.up2, 'scaled:scale') + attr(xi.models$y.train.up2, 'scaled:center')
    xi.dose.predict.up[d1.d2.ind$data2.ind] <- c(xi.dose.predict.up2) * attr(xi.models$y.train.up1, 'scaled:scale') + attr(xi.models$y.train.up1, 'scaled:center')
    
  } else {
    stop('xi.method must be "regression_forest" or "neural_network"')
  }
  
  
  # For each bandwidth/window
  PEI.bounds <- foreach(i=1:windows.length, .combine="rbind") %do% {
    
    # Compute the kernel estimations
    kernel.est.down <- K.h(dataset$t - dose, windows[i])
    kernel.est.up <- kernel.est.down  # Keep same kernel estimation for lower and upper bounds
    
    data.min.idx <- dataset$Y <= Q.predict[, 1]  # q_{1-gamma}^{X_i, T_i}
    data.max.idx <- dataset$Y > Q.predict[, 2]  # q_{gamma}^{X_i, T_i}
    
    if (stabilization) {
      stab.min <- mean((kernel.est.down[data.min.idx] / cond.density[data.min.idx]))
      stab.min.unconf <- mean((kernel.est.down[data.min.idx] / cond.density[data.min.idx]))
      stab.max <- mean((kernel.est.up[data.max.idx] / cond.density[data.max.idx]))
      stab.max.unconf <- mean((kernel.est.up[data.max.idx] / cond.density[data.max.idx]))
    } else {
      stab.max <- 1
      stab.max.unconf <- 1
      stab.min <- 1
      stab.min.unconf <- 1
    }
    
    tau <- lambda / (1 + lambda)
    
    if (is.null(xi.method)) {  # If no method for xi was provided, compute simple estimators
      
      # Minimization problem (lower bound)
      minimum <- mean(Y.predict.dose) + ((2*tau-1)/tau) * mean((kernel.est.down[data.min.idx] * (dataset$Y[data.min.idx] - Y.predict.all[data.min.idx]) / cond.density[data.min.idx])) / stab.min
      minimum.unconf <- mean(Y.predict.dose) + ((2*tau-1)/tau) * mean((kernel.est.down[data.min.idx] * (dataset$Y[data.min.idx] - Y.predict.all[data.min.idx]) / cond.density[data.min.idx])) / stab.min.unconf
      
      # Maximization problem (upper bound)
      maximum <- mean(Y.predict.dose) + ((2*tau-1)/tau) * mean((kernel.est.up[data.max.idx] * (dataset$Y[data.max.idx] - Y.predict.all[data.max.idx]) / cond.density[data.max.idx])) / stab.max
      maximum.unconf <- mean(Y.predict.dose) + ((2*tau-1)/tau) * mean((kernel.est.up[data.max.idx] * (dataset$Y[data.max.idx] - Y.predict.all[data.max.idx]) / cond.density[data.max.idx])) / stab.max.unconf
      
    } else { # Else, compute doubly robust estimators
      
      minimum <- mean(xi.dose.predict.down) + mean(kernel.est.down * (dataset$Y * lambda**(-sign(dataset$Y - Q.predict[, 1])) - xi.t.predict.down) / cond.density) / mean(kernel.est.down / cond.density)  # Essayer de stabiliser
      minimum.unconf <- minimum  # Finally, this is not used
      maximum <- mean(xi.dose.predict.up) + mean(kernel.est.up * (dataset$Y * lambda**(sign(dataset$Y - Q.predict[, 2])) - xi.t.predict.up) / cond.density) / mean(kernel.est.up / cond.density)
      maximum.unconf <- maximum  # Finally, this is not used
      
    }
    
    c(minimum, maximum, minimum.unconf, maximum.unconf)
  }
  
  return(PEI.bounds)
}


# The neural network used to compute the generalized propensity score (GPS) p(t|x)
base_neural_network_gps <- nn_module(classname="BaseNN_GPS",
                                     initialize=function(dim_cov, dim_output, num_components, dim_hidden) {
                                       
                                       self$feature_extractor <- nn_sequential(
                                         nn_linear(dim_cov, dim_hidden, bias=T),
                                         nn_leaky_relu(negative_slope=0.04),
                                         nn_dropout(p=0.04),
                                         nn_linear(dim_hidden, dim_hidden, bias=T),
                                         nn_leaky_relu(negative_slope=0.04),
                                         nn_dropout(p=0.04)
                                       )
                                       
                                       self$density_estimator <- nn_sequential(
                                         nn_linear(dim_hidden, dim_hidden*2, bias=T),
                                         nn_leaky_relu(negative_slope=0.04),
                                         nn_dropout(p=0.04),
                                         nn_linear(dim_hidden*2, dim_hidden*2, bias=T),
                                         nn_leaky_relu(negative_slope=0.04),
                                         nn_dropout(p=0.04),
                                         gaussian_mixture_regression(dim_hidden*2, dim_output, num_components))
                                     },
                                     forward=function(covariates) {
                                       phi <- self$feature_extractor(covariates)  # Of dim n.obs*dim_hidden
                                       return(self$density_estimator(phi))
                                     })

