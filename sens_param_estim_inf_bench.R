# Import libraries
library(foreach)
library(torch)
library(latex2exp)  # For LaTeX expressions
library(tictoc)
library(ggplot2)

# Import R files
source("./utils.R")


data.name.arg <- "cmr"
version.arg <- 1


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
  n.cov <- ncol(X)
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

max.iter.gps <- 1000
patience.gps <- 5

fine.tun.nn.params <- list(train.prop.gps=0.8, valid.prop.gps=0.1, n.random.splits.gps=2, max.iter.gps=max.iter.gps, patience.gps=patience.gps)
nn.params <- list(max.iter.gps=max.iter.gps, patience.gps=patience.gps)
nn.init <- NULL

# Function that estimates the GPS with 2-fold cross-fitting
gps.fun <- function(X, Y, t, data.name,
                    nn.init=NULL, D1.prop=0.6, fine.tun.nn.params=list(train.prop.gps=0.8, valid.prop.gps=0.1, n.random.splits.gps=2, max.iter.gps=1000, patience.gps=20),
                    nn.params=list(max.iter.gps=1000, patience.gps=20),
                    grid.K=NULL, grid.hid.dim=NULL, grid.lr=NULL,
                    device=torch_device("cpu"), verbose=TRUE) {
  
  n.all <- nrow(X)
  
  # Be careful here with the dimensions of the tensors!
  X.tensor <- torch_tensor(as.matrix(X), device=device)
  t.tensor <- torch_tensor(c(t), device=device)
  y.tensor <- torch_tensor(c(Y), device=device)
  
  # Divide dataset D into D1 and D2
  n.data1 <- floor(D1.prop*n.all)
  n.data2 <- n.all - n.data1
  data1.ind <- sort(sample(1:n.all, n.data1, replace=FALSE))
  data2.ind <- setdiff(1:n.all, data1.ind)  # Because X.tensor[-data1.ind, ] does not work
  # Store the indices of D1 and D2
  d1.d2.ind <- list(data1.ind=data1.ind, data2.ind=data2.ind)
  
  X1 <- X[data1.ind, ]
  X2 <- X[data2.ind, ]
  t1 <- t[data1.ind]
  t2 <- t[data2.ind]
  Y1 <- Y[data1.ind]
  Y2 <- Y[data2.ind]
  
  # Get corresponding tensors
  X.tensor1 <- X.tensor[data1.ind, ]
  t.tensor1 <- t.tensor[data1.ind]
  y.tensor1 <- y.tensor[data1.ind]
  
  X.tensor2 <- X.tensor[data2.ind, ]
  t.tensor2 <- t.tensor[data2.ind]
  y.tensor2 <- y.tensor[data2.ind]
  
  # Divide D1 into train and validation sets
  train.valid.data1 <- train.valid.split(n.data=n.data1, X.tensor=X.tensor1,
                                         t.tensor=t.tensor1, y.tensor=y.tensor1)
  # Divide D2 into train and validation sets
  train.valid.data2 <- train.valid.split(n.data=n.data2, X.tensor=X.tensor2,
                                         t.tensor=t.tensor2, y.tensor=y.tensor2)
  
  if (is.null(nn.init)) {
    
    # If no model was given, fine-tune and train neural networks from scratch on D1
    nn.init <- list(gps.gmm=NULL, resp.gmm=NULL)
    
    gps.file.name <- paste0("./params/", data.name, "/optimal_params_gps.RData")
    
    if (file.exists(gps.file.name)) {
      
      # Get optimal parameters if they are already stored
      optimal.params.gps <- readRDS(gps.file.name)
      
    } else {
      
      # Fine-tuning for the Generalized Propensity Score p(t|X)
      # It is ok to get the fine-tuned parameters only on D1
      optimal.params.gps <- nn.fine.tuning(K.vec=grid.K,
                                           dim.hidden.vec=grid.hid.dim,
                                           lr.vec=grid.lr,
                                           X.tensor=X.tensor1,
                                           t.tensor=t.tensor1,
                                           train.prop=fine.tun.nn.params$train.prop.gps,
                                           valid.prop=fine.tun.nn.params$valid.prop.gps,
                                           n.random.splits=fine.tun.nn.params$n.random.splits.gps,
                                           max.iter=fine.tun.nn.params$max.iter.gps,
                                           patience=fine.tun.nn.params$patience.gps,
                                           device=device,
                                           verbose=verbose)
      
      saveRDS(optimal.params.gps, file=gps.file.name)
      
    }
    
    # Train the final model on 90% of D1 and validate on 10% of D1
    trained.gps.model1 <- train.nn(nn.architecture=base_neural_network_gps,
                                   X.tensor.train=train.valid.data1$X.tensor.train,
                                   X.tensor.valid=train.valid.data1$X.tensor.valid,
                                   t.tensor.train=train.valid.data1$t.tensor.train,
                                   t.tensor.valid=train.valid.data1$t.tensor.valid,
                                   max.iter=nn.params$max.iter.gps,
                                   patience=nn.params$patience.gps,
                                   K=optimal.params.gps$K.optim,
                                   lr=optimal.params.gps$lr.optim,
                                   dim.hidden=optimal.params.gps$dim.hidden.optim,
                                   device=device,
                                   verbose=verbose)
    
    nn.init$gps.gmm1 <- trained.gps.model1$gmm
    
    # Train the final model on 90% of D2 and validate on 10% of D2
    trained.gps.model2 <- train.nn(nn.architecture=base_neural_network_gps,
                                   X.tensor.train=train.valid.data2$X.tensor.train,
                                   X.tensor.valid=train.valid.data2$X.tensor.valid,
                                   t.tensor.train=train.valid.data2$t.tensor.train,
                                   t.tensor.valid=train.valid.data2$t.tensor.valid,
                                   max.iter=nn.params$max.iter.gps,
                                   patience=nn.params$patience.gps,
                                   K=optimal.params.gps$K.optim,
                                   lr=optimal.params.gps$lr.optim,
                                   dim.hidden=optimal.params.gps$dim.hidden.optim,
                                   device=device,
                                   verbose=verbose)
    
    nn.init$gps.gmm2 <- trained.gps.model2$gmm
    
    
  }
  
  # Put the networks in evaluation mode
  nn.init$gps.gmm1$eval()
  nn.init$gps.gmm2$eval()
  
  # Evaluate the densities on D1 and D2
  gps.eval1 <- as.array(exp(nn.init$gps.gmm2(covariates=X.tensor1)$log_prob(x=t.tensor1)$to(device="cpu")))
  gps.eval2 <- as.array(exp(nn.init$gps.gmm1(covariates=X.tensor2)$log_prob(x=t.tensor2)$to(device="cpu")))
  
  # Store the evaluations
  gps.eval.all <- rep(NA, n.all)
  gps.eval.all[data1.ind] <- gps.eval1
  gps.eval.all[data2.ind] <- gps.eval2
  
  return(list(gps.eval.all=gps.eval.all,
              X2=X2, t2=t2, Y2=Y2,
              d1.d2.ind=d1.d2.ind,
              nn.init=nn.init,
              optimal.params.gps=optimal.params.gps))
}


# Vector that contains the estimated gammas for each observed confounder
est.gammas.vec <- rep(NA, n.cov)
est.gammas <- NULL
clipping.ind <- NULL
method.vec <- NULL

# Compute the GPS conditionally on all observed covariates. This will be the numerator.
gps.list.all.conf <- gps.fun(X=X, Y=Y, t=t, data.name=data.name, nn.init=NULL,
                             D1.prop=D1.prop, fine.tun.nn.params=fine.tun.nn.params,
                             nn.params=nn.params, grid.K=K.vec, grid.hid.dim=dim.hidden.vec, grid.lr=rd.lr.ind,
                             device=device, verbose=verbose)

gps.all.conf <- gps.list.all.conf$gps.eval.all

# Remove one-by-one each observed covariate and compute the GPS conditionally on all covariates except the one that was removed. This will be the denominator.
for (i in 1:n.cov) {
  message(paste0("-- Confounder ", i, ": ", cov.names[i], " --"))

  # Remove covariate i from the observed covariates
  X.minus.i <- X[, -i]

  gps.list.conf.minus.i <- gps.fun(X=X.minus.i, Y=Y, t=t, data.name=data.name, nn.init=NULL,
                                   D1.prop=D1.prop, fine.tun.nn.params=fine.tun.nn.params,
                                   nn.params=nn.params, grid.K=K.vec, grid.hid.dim=dim.hidden.vec, grid.lr=rd.lr.ind,
                                   device=device, verbose=verbose)
  
  gps.conf.minus.i <- gps.list.conf.minus.i$gps.eval.all
  
  message("* Without weight clipping *")
  
  estim.gammas <- gamma.estim.fun(gps.all.conf=gps.all.conf,
                                  gps.conf.minus.i=gps.conf.minus.i)
  
  message("* With weight clipping *")
  
  # If we do propensity score clipping
  gps.all.conf.clipped <- gps.all.conf
  gps.all.conf.clipped[gps.all.conf.clipped < 0.1] <- 0.1
  gps.conf.minus.i.clipped <- gps.conf.minus.i
  gps.conf.minus.i.clipped[gps.conf.minus.i.clipped < 0.1] <- 0.1
  
  estim.gammas.clipped <- gamma.estim.fun(gps.all.conf=gps.all.conf.clipped,
                                          gps.conf.minus.i=gps.conf.minus.i.clipped)
  
  est.gammas.vec[i] <- estim.gammas.clipped$estimated.gamma
}
