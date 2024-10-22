# All functions used by the algorithm from Jesson et al. and us

# Import functions used particularly by each method
source("jesson_fun.R")
source("cont_qb_fun.R")
source("simulated_data_fun.R")


# Function to preprocess the data
# Returns a list with data.frame object (all.data), number of individuals, covariate names, scaled outcome and scaled treatment
preprocess.pm2.5.cmr.data <- function(folder.name) {
  
  # Import PM2.5 and Cardiovascular Mortality Rate (CMR) data
  pm25.cmr.data <- read.csv(paste(folder.name, "County_annual_PM25_CMR.csv", sep=""))[, -1]  # Remove first column which is only an index
  raw.variables.data <- read.csv(paste(folder.name, "County_RAW_variables.csv", sep=""))[, -1]
  ses.index.quintile.data <- read.csv(paste(folder.name, "County_SES_index_quintile.csv", sep=""))[, -1]
  
  # Select only relevant columns
  raw.variables.col.interest <- c("FIPS", "healthfac_2005_1999", "population_2000",
                                  "civil_unemploy_2010", "median_HH_inc_2010",
                                  "femaleHH_ns_pct_2010", "vacant_HHunit_2010",
                                  "owner_occ_pct_2010", "eduattain_HS_2010",
                                  "pctfam_pover_2010")
  ses.index.col.interest <- c("FIPS", "SES_index_2010")
  
  # Get covariate names
  cov.names <- c(raw.variables.col.interest[-1], ses.index.col.interest[-1])
  
  # Get covariates
  raw.variables.data <- raw.variables.data[, raw.variables.col.interest]
  ses.index.quintile.data <- ses.index.quintile.data[, ses.index.col.interest]
  
  # Keep only year 2010
  pm25.cmr.data <- pm25.cmr.data[pm25.cmr.data$Year == 2010, ]
  pm25.cmr.data <- subset(pm25.cmr.data, select=-c(Year))  # Remove Year column
  
  # Merge dataframes
  all.data <- merge(pm25.cmr.data, raw.variables.data, by="FIPS", all.x=TRUE)
  all.data <- merge(all.data, ses.index.quintile.data, by="FIPS", all.x=TRUE)
  
  # Get number of data
  n.all <- nrow(all.data)
  
  # Save unnormalized data
  unnormalized.data <- all.data
  
  # Log-normalize population_2000 and median_HH_inc_2010
  all.data$population_2000 <- log(all.data$population_2000)
  all.data$median_HH_inc_2010 <- log(all.data$median_HH_inc_2010)
  
  # Center and scale all covariates
  all.data[, cov.names] <- scale(all.data[, cov.names])
  
  # Center and scale CMR and PM2.5
  scaled.Y <- scale(all.data$CMR)  # Here, we can save the mean and standard deviation
  scaled.t <- scale(all.data$PM2.5)  # Here, we can save the mean and standard deviation
  
  all.data$CMR <- scaled.Y[, 1]
  all.data$PM2.5 <- scaled.t[, 1]
  
  # Remove 10% of the outliers as in the simulated data
  hat.values <- hatvalues(lm(rnorm(n.all, 0, 1) ~ all.data$PM2.5 + as.matrix(all.data[, cov.names]) + all.data$CMR))
  outliers <- which(hat.values > quantile(hat.values, 0.9))
  n.all <- n.all - length(outliers)
  
  # Rename the columns of the outcome Y and the treatment t
  names(all.data)[names(all.data) == "CMR"] <- "Y"
  names(all.data)[names(all.data) == "PM2.5"] <- "t"
  
  # Remove column FIPS
  all.data <- subset(all.data, select=-c(FIPS))
  
  # Remove the outliers
  all.data <- all.data[-outliers, ]
  unnormalized.data <- unnormalized.data[-outliers, ]
  
  return(list(unnormalized.data=unnormalized.data, normalized.data=all.data, n.all=n.all, cov.names=cov.names, scaled.Y=scaled.Y, scaled.t=scaled.t))
}


# Function to delete more easily some files
delete.optimal.params.files <- function(data.name) {
  
  gps.file.name <- paste0("./params/", data.name, "/optimal_params_gps.RData")
  resp.file.name <- paste0("./params/", data.name, "/optimal_params_resp.RData")
  
  if (file.exists(gps.file.name)) {
    file.remove(gps.file.name)
    message(paste(gps.file.name, "successfully removed"))
  } else {
    message(paste(gps.file.name, "does not exist"))
  }
  
  if (file.exists(resp.file.name)) {
    file.remove(resp.file.name)
    message(paste(resp.file.name, "successfully removed"))
  } else {
    message(paste(resp.file.name, "does not exist"))
  }
  
}


# The neural network used to compute p(y|x,t)
base_neural_network <- nn_module(classname="BaseNN",
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
                                     nn_linear(dim_hidden+1, dim_hidden*2, bias=T),  # Be careful, dim_cov+1 for the treatment
                                     nn_leaky_relu(negative_slope=0.04),
                                     nn_dropout(p=0.04),
                                     nn_linear(dim_hidden*2, dim_hidden*2, bias=T),
                                     nn_leaky_relu(negative_slope=0.04),
                                     nn_dropout(p=0.04),
                                     gaussian_mixture_regression(dim_hidden*2, dim_output, num_components))
                                 },
                                 forward=function(covariates, treatment) {
                                   treatment <- torch_unsqueeze(treatment, -1)  # Of dim n.obs*1
                                   phi <- self$feature_extractor(covariates)  # Of dim n.obs*dim_hidden
                                   phi_cat <- torch_cat(c(phi, treatment), dim=-1)  # Of dim n.obs*(dim_hidden+1)
                                   return(self$density_estimator(phi_cat))
                                 })


# Create a module of name GMR with function initialize and forward
gaussian_mixture_regression <- nn_module(classname="GMR",
                                         
                                         # num_components is the number of components we want in the GMM
                                         initialize=function(dim_input, dim_output, num_components) {
                                           num_out_mu_sigma <- num_components * dim_output
                                           
                                           self$mu <- nn_linear(
                                             in_features=dim_input, out_features=num_out_mu_sigma, bias=T
                                           )  # We want one mean for each component
                                           
                                           sigma <- nn_linear(
                                             in_features=dim_input, out_features=num_out_mu_sigma, bias=T
                                           )  # We want one standard deviation for each component
                                           
                                           self$pi <- nn_linear(
                                             in_features=dim_input, out_features=num_components, bias=T
                                           )  # We want one weight for each component
                                           
                                           self$sigma <- nn_sequential(sigma, nn_softplus())  # Softplus activation to have positive values
                                           self$num_components <- num_components
                                           self$dim_output <- dim_output
                                         },
                                         
                                         forward=function(inputs) {
                                           loc <- self$mu(inputs)
                                           scale <- self$sigma(inputs) + 1e-6
                                           mixture_distribution <- distr_categorical(logits=self$pi(inputs))
                                           component_distribution <- distr_normal(loc=loc, scale=scale)
                                           # distr_mixture_same_family performs everything to sample from the estimated distribution and allows backpropagation even if there is randomness
                                           gmm_distr <- distr_mixture_same_family(mixture_distribution, component_distribution)
                                           
                                           # We want to minimize the negative log-likelihood
                                           return(gmm_distr)
                                         }
)


train.nn <- function(nn.architecture, nn.model=NULL,
                     X.tensor.train, X.tensor.valid, X.tensor.test=NULL,
                     t.tensor.train, t.tensor.valid, t.tensor.test=NULL,
                     y.tensor.train=NULL, y.tensor.valid=NULL, y.tensor.test=NULL,
                     max.iter=1000, patience=40, patience.range=5,
                     K=3, lr=1e-3, dim.hidden=16, device, verbose=TRUE) {
  
  # Initialize the model
  if (!is.null(y.tensor.train)) {  # For p(y|x,t)

    if (is.null(nn.model)) {
      gmm <- nn.architecture(dim_cov=X.tensor.train$size(2), dim_output=1, num_components=K, dim_hidden=dim.hidden)$to(device=device)  # Arguments passed to initialize
    } else {
      # If a model was given, initialize a new model and load the weights
      # Transfer learning
      gmm <- nn.architecture(dim_cov=X.tensor.train$size(2), dim_output=1, num_components=K, dim_hidden=dim.hidden)$to(device=device)
      gmm$load_state_dict(nn.model$state_dict(prefix=""))
    }
    
  } else {  # For the GPS p(t|x)

    if (is.null(nn.model)) {
      gmm <- nn.architecture(dim_cov=X.tensor.train$size(2), dim_output=1, num_components=K, dim_hidden=dim.hidden)$to(device=device)  # Arguments passed to initialize
    } else {
      # If a model was given
      gmm <- nn.architecture(dim_cov=X.tensor.train$size(2), dim_output=1, num_components=K, dim_hidden=dim.hidden)$to(device=device)
      gmm$load_state_dict(nn.model$state_dict(prefix=""))
    }
    
  }

  # Initialize the optimizer (Adam)
  optim <- optim_adam(gmm$parameters, lr=lr)
  
  # Initialize the training loss and validation loss curves
  learning.curve.train <- NULL
  learning.curve.valid <- NULL
  
  # Initialize vectors to save patience.range validation losses at time t and t+patience
  saved.valid.losses <- rep(NA, patience.range)
  patience.valid.losses <- rep(NA, patience.range)
  
  if (verbose) {
    # Initialize the progress bar
    pb <- txtProgressBar(min=0,      # Minimum value of the progress bar
                         max=max.iter, # Maximum value of the progress bar
                         style=3,    # Progress bar style (also available style = 1 and style = 2)
                         width=50,   # Progress bar width. Defaults to getOption("width")
                         char="=")   # Character used to create the bar
  }
  
  for (i in 0:(max.iter-1)) {
    optim$zero_grad()
    
    if (!is.null(y.tensor.train)) {  # For p(y|x,t)
      loss.train <- -gmm(covariates=X.tensor.train, treatment=t.tensor.train)$log_prob(y.tensor.train)$mean()  # Training loss
      loss.valid <- -gmm(covariates=X.tensor.valid, treatment=t.tensor.valid)$log_prob(y.tensor.valid)$mean()  # Validation loss
    } else {  # For the GPS p(t|x)
      loss.train <- -gmm(covariates=X.tensor.train)$log_prob(t.tensor.train)$mean()  # Training loss
      loss.valid <- -gmm(covariates=X.tensor.valid)$log_prob(t.tensor.valid)$mean()  # Validation loss
    }
    
    saved.valid.losses[i%%patience.range+1] <- as.numeric(loss.valid$to(device="cpu"))
    
    # Initialize mean.saved.valid.loss
    if (i == patience.range-1) {
      mean.saved.valid.loss <- mean(saved.valid.losses)
    }
    
    loss.train$backward()
    optim$step()
    learning.curve.train <- c(learning.curve.train, loss.train$item())
    learning.curve.valid <- c(learning.curve.valid, loss.valid$item())
    
    if (verbose) {
      setTxtProgressBar(pb, i+1)  # Add one unit to the progress bar
    }
    
    if ((i%%patience+1 == patience) & (i%%patience.range+1 == patience.range)) {
      
      mean.current.valid.loss <- mean(saved.valid.losses)
      
      # Stop the loop if, after patience iterations, the mean of the validation losses increased
      if (mean.current.valid.loss > mean.saved.valid.loss) {
        if (verbose) {
          close(pb)
        }
        break
      }
      
      mean.saved.valid.loss <- mean.current.valid.loss
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  if (!is.null(y.tensor.train)) {  # For p(y|x,t)
    
    # If both X.tensor.test, t.tensor.test and y.tensor.test are not null
    if (!is.null(X.tensor.test) & !is.null(t.tensor.test) & !is.null(y.tensor.test)) {
      # Compute the loss on the test set
      test.loss <- as.numeric(-gmm(covariates=X.tensor.test, treatment=t.tensor.test)$log_prob(y.tensor.test)$mean())
    } else {
      test.loss <- NULL
    }
    
  } else {  # For the GPS p(t|x)
    
    # If both X.tensor.test and t.tensor.test are not null
    if (!is.null(X.tensor.test) & !is.null(t.tensor.test)) {
      # Compute the loss on the test set
      test.loss <- as.numeric(-gmm(covariates=X.tensor.test)$log_prob(t.tensor.test)$mean())
    } else {
      test.loss <- NULL
    }
    
  }

  return(list(gmm=gmm, test.loss=test.loss, learning.curve.train=learning.curve.train, learning.curve.valid=learning.curve.valid))
}


# Function to fine-tune the hyperparameters of the neural networks
nn.fine.tuning <- function(K.vec, dim.hidden.vec, lr.vec,
                           X.tensor, t.tensor, y.tensor=NULL,
                           train.prop=0.8, valid.prop=0.1,
                           n.random.splits=2, max.iter=1500,
                           patience=40, patience.range=5,
                           device,
                           verbose=TRUE) {
  
  # Length of search space
  search.space.len <- length(lr.vec)
  # Number of data
  n.data <- nrow(X.tensor)
  
  # Get train, validation and test sets sample size
  n.train <- floor(train.prop*n.data)
  n.valid <- floor(valid.prop*n.data)
  n.test <- n.data - n.train - n.valid
  
  # Initialize the matrix that will contain the test losses
  test.loss.mat <- matrix(nrow=search.space.len, ncol=n.random.splits)
  # Initialize the lists that will contain the loss curves
  lc.train.list <- list()
  lc.valid.list <- list()
  
  for (i in 1:n.random.splits) {
    
    # Initialize the lists that will contain the loss curves
    lc.train.list.i <- list()
    lc.valid.list.i <- list()
    
    # Random split of D1 into train, validation and test
    train.ind <- sample(1:n.data, n.train)
    valid.ind <- sample(setdiff(1:n.data, train.ind), n.valid)
    test.ind <- setdiff(1:n.data, c(train.ind, valid.ind))
    
    X.tensor.train <- X.tensor[train.ind, ]
    X.tensor.valid <- X.tensor[valid.ind, ]
    X.tensor.test <- X.tensor[test.ind, ]
    
    t.tensor.train <- t.tensor[train.ind]
    t.tensor.valid <- t.tensor[valid.ind]
    t.tensor.test <- t.tensor[test.ind]
    
    # If y.tensor is given, split also into train, validation and test
    if (!is.null(y.tensor)) {
      y.tensor.train <- y.tensor[train.ind]
      y.tensor.valid <- y.tensor[valid.ind]
      y.tensor.test <- y.tensor[test.ind]
    }
    
    for (j in 1:search.space.len) {
      
      # If y.tensor was not given, train the GPS (p(t|x)) neural network
      if (is.null(y.tensor)) {
        
        trained.model <- train.nn(nn.architecture=base_neural_network_gps,
                                  X.tensor.train=X.tensor.train,
                                  X.tensor.valid=X.tensor.valid,
                                  X.tensor.test=X.tensor.test,
                                  t.tensor.train=t.tensor.train,
                                  t.tensor.valid=t.tensor.valid,
                                  t.tensor.test=t.tensor.test,
                                  max.iter=max.iter, patience=patience,
                                  patience.range=patience.range,
                                  K=K.vec[j], lr=lr.vec[j],
                                  dim.hidden=dim.hidden.vec[j],
                                  device=device,
                                  verbose=verbose)
        
      } else {  # Else, train the NN to estimate p(y|x,t)
        
        trained.model <- train.nn(nn.architecture=base_neural_network,
                                  X.tensor.train=X.tensor.train,
                                  X.tensor.valid=X.tensor.valid,
                                  X.tensor.test=X.tensor.test,
                                  t.tensor.train=t.tensor.train,
                                  t.tensor.valid=t.tensor.valid,
                                  t.tensor.test=t.tensor.test,
                                  y.tensor.train=y.tensor.train,
                                  y.tensor.valid=y.tensor.valid,
                                  y.tensor.test=y.tensor.test,
                                  max.iter=max.iter, patience=patience,
                                  patience.range=patience.range,
                                  K=K.vec[j], lr=lr.vec[j],
                                  dim.hidden=dim.hidden.vec[j],
                                  device=device,
                                  verbose=verbose)
        
      }
      
      test.loss.mat[j, i] <- trained.model$test.loss
      lc.train.list.i[[j]] <- trained.model$learning.curve.train
      lc.valid.list.i[[j]] <- trained.model$learning.curve.valid
    }
    
    lc.train.list[[i]] <- lc.train.list.i
    lc.valid.list[[i]] <- lc.valid.list.i
  }
  
  # Compute the mean test loss for each hyperparameter combination
  mean.loss.vec <- rowMeans(test.loss.mat)
  # Get the index of the minimum loss
  optimal.ind <- which.min(mean.loss.vec)
  
  # Get the optimal hyperparameters
  K.optim <- K.vec[optimal.ind]
  dim.hidden.optim <- dim.hidden.vec[optimal.ind]
  lr.optim <- lr.vec[optimal.ind]
  
  return(list(K.optim=K.optim, dim.hidden.optim=dim.hidden.optim,
              lr.optim=lr.optim, mean.loss.vec=mean.loss.vec,
              lc.train.list=lc.train.list, lc.valid.list=lc.valid.list))
}


# Function to compute the PEI on a sample (X, Y, t)
PEI.fun <- function(X, Y, t, data.name,
                    doses, gamma, bootstrap.ind=NULL, stabilization=TRUE,
                    bandwidths, B.param=50,
                    compute.QB=TRUE, cond.quant.method=c("quantile_forest"), Q.predict=NULL, xi.method=c("neural_network"),
                    compute.Jesson=TRUE, X.sample.len=NULL, Y.sample.len=2000,
                    nn.init=NULL, xi.models=NULL, D1.prop=0.6, fine.tun.nn.params=list(train.prop.gps=0.8, valid.prop.gps=0.1, n.random.splits.gps=2, max.iter.gps=1000, patience.gps=20),
                    nn.params=list(max.iter.gps=1000, max.iter.resp=1000, patience.gps=20, patience.resp=20),
                    grid.K=NULL, grid.hid.dim=NULL, grid.lr=NULL, use.parallel.jesson=TRUE, device=torch_device("cpu"), verbose=TRUE) {
  
  start.time.qb1 <- Sys.time()
  start.time.jesson1 <- Sys.time()
  
  tau <- gamma / (1 + gamma)
  doses.length <- length(doses)
  
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
  
  end.time.qb1 <- Sys.time()
  end.time.jesson1 <- Sys.time()
  
  if (is.null(nn.init)) {
    
    start.time.qb2 <- Sys.time()
    start.time.jesson2 <- Sys.time()
    
    # If no model was given, fine-tune and train neural networks from scratch on D1
    nn.init <- list(gps.gmm=NULL, resp.gmm=NULL)
    
    resp.file.name <- paste0("./params/", data.name, "/optimal_params_resp.RData")
    
    if (file.exists(resp.file.name)) {
      
      # Get optimal parameters if they are already stored
      optimal.params.resp <- readRDS(resp.file.name)
      
    } else {
      
      # Fine-tuning for the response density p(Y|X,t)
      optimal.params.resp <- nn.fine.tuning(K.vec=grid.K,
                                            dim.hidden.vec=grid.hid.dim,
                                            lr.vec=grid.lr,
                                            X.tensor=X.tensor1,
                                            t.tensor=t.tensor1,
                                            y.tensor=y.tensor1,
                                            train.prop=fine.tun.nn.params$train.prop.resp,
                                            valid.prop=fine.tun.nn.params$valid.prop.resp,
                                            n.random.splits=fine.tun.nn.params$n.random.splits.resp,
                                            max.iter=fine.tun.nn.params$max.iter.resp,
                                            patience=fine.tun.nn.params$patience.resp,
                                            device=device,
                                            verbose=verbose)
      
      saveRDS(optimal.params.resp, file=resp.file.name)
      
    }
    
    # Train the final model on 90% of D1 and validate on 10% of D1
    trained.resp.model1 <- train.nn(nn.architecture=base_neural_network,
                                    X.tensor.train=train.valid.data1$X.tensor.train,
                                    X.tensor.valid=train.valid.data1$X.tensor.valid,
                                    t.tensor.train=train.valid.data1$t.tensor.train,
                                    t.tensor.valid=train.valid.data1$t.tensor.valid,
                                    y.tensor.train=train.valid.data1$y.tensor.train,
                                    y.tensor.valid=train.valid.data1$y.tensor.valid,
                                    max.iter=nn.params$max.iter.resp,
                                    patience=nn.params$patience.resp,
                                    K=optimal.params.resp$K.optim,
                                    lr=optimal.params.resp$lr.optim,
                                    dim.hidden=optimal.params.resp$dim.hidden.optim,
                                    device=device,
                                    verbose=verbose)
    
    nn.init$resp.gmm1 <- trained.resp.model1$gmm
    
    # Train the final model on 90% of D1 and validate on 10% of D1
    trained.resp.model2 <- train.nn(nn.architecture=base_neural_network,
                                    X.tensor.train=train.valid.data2$X.tensor.train,
                                    X.tensor.valid=train.valid.data2$X.tensor.valid,
                                    t.tensor.train=train.valid.data2$t.tensor.train,
                                    t.tensor.valid=train.valid.data2$t.tensor.valid,
                                    y.tensor.train=train.valid.data2$y.tensor.train,
                                    y.tensor.valid=train.valid.data2$y.tensor.valid,
                                    max.iter=nn.params$max.iter.resp,
                                    patience=nn.params$patience.resp,
                                    K=optimal.params.resp$K.optim,
                                    lr=optimal.params.resp$lr.optim,
                                    dim.hidden=optimal.params.resp$dim.hidden.optim,
                                    device=device,
                                    verbose=verbose)
    
    nn.init$resp.gmm2 <- trained.resp.model2$gmm
    
    end.time.jesson2 <- Sys.time()
    
    
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
    
    if (compute.QB) {
      
      if (cond.quant.method == "quantile_forest") {
        
        Q.model <- quantile_forest(X=data.frame(X1, t=t1), Y=Y1,
                                   quantiles=c(1-tau, tau))
        
        # Estimate the conditional quantiles Q_tau(Y|X_i, dose) and Q_{1-tau}(Y|X, dose)
        Q.predict <- predict(Q.model, data.frame(X, t),
                             quantiles=c(1-tau, tau))$predictions
        
      } else if (cond.quant.method == "root_search") {
        
        # Estimate the conditional quantiles using the fitted distribution of p(y|x,t)
        # On D1, using the network fitted on D2
        nn.init$Q.predict1 <- cond.quant.estim(X=X.tensor1, t=t.tensor1,
                                               resp.gmm=nn.init$resp.gmm2,
                                               tau=tau, device=device, verbose=verbose)
        # On D2, using the network fitted on D1
        nn.init$Q.predict2 <- cond.quant.estim(X=X.tensor2, t=t.tensor2,
                                               resp.gmm=nn.init$resp.gmm1,
                                               tau=tau, device=device, verbose=verbose)
        
        # This is stored for the bootstrap resamples
        Q.predict <- matrix(NA, nrow=n.all, ncol=2)
        Q.predict[data1.ind, ] <- nn.init$Q.predict1
        Q.predict[data2.ind, ] <- nn.init$Q.predict2
        
      } else {
        stop("cond.quant.method must be 'quantile_forest' or 'root_search'")
      }
      
      
      if (is.null(xi.method)) {
        
        xi.models <- NULL
        
      } else if (xi.method == "regression_forest") {
        
        # Train function xi that brings double-robustness
        # Model fitting on D1
        X.t.d1 <- cbind(X1, t1)
        names(X.t.d1)[names(X.t.d1) == "t1"] <- "t"
        X.t.d1 <- as.matrix(X.t.d1)
        # For the lower bound
        Y.gamma.down.d1 <- scale(Y1 * gamma**(-sign(Y1 - nn.init$Q.predict1[, 1])))
        xi.model.down.d1 <- regression_forest(X=X.t.d1,
                                              Y=Y.gamma.down.d1,
                                              num.trees=100,
                                              tune.parameters="all")
        # For the upper bound
        Y.gamma.up.d1 <- scale(Y1 * gamma**(sign(Y1 - nn.init$Q.predict1[, 2])))
        xi.model.up.d1 <- regression_forest(X=X.t.d1,
                                            Y=Y.gamma.up.d1,
                                            num.trees=100,
                                            tune.parameters="all")
        
        # Model fitting on D2
        X.t.d2 <- cbind(X2, t2)
        names(X.t.d2)[names(X.t.d2) == "t2"] <- "t"
        X.t.d2 <- as.matrix(X.t.d2)
        # For the lower bound
        Y.gamma.down.d2 <- scale(Y2 * gamma**(-sign(Y2 - nn.init$Q.predict2[, 1])))
        xi.model.down.d2 <- regression_forest(X=X.t.d2,
                                              Y=Y.gamma.down.d2,
                                              num.trees=100,
                                              tune.parameters="all")
        # For the upper bound
        Y.gamma.up.d2 <- scale(Y2 * gamma**(sign(Y2 - nn.init$Q.predict2[, 2])))
        xi.model.up.d2 <- regression_forest(X=X.t.d2,
                                            Y=Y.gamma.up.d2,
                                            num.trees=100,
                                            tune.parameters="all")
        
        xi.models <- list(xi.model.down.d1=xi.model.down.d1,
                          xi.model.up.d1=xi.model.up.d1,
                          xi.model.down.d2=xi.model.down.d2,
                          xi.model.up.d2=xi.model.up.d2,
                          Y.gamma.down.d1=Y.gamma.down.d1,
                          Y.gamma.up.d1=Y.gamma.up.d1,
                          Y.gamma.down.d2=Y.gamma.down.d2,
                          Y.gamma.up.d2=Y.gamma.up.d2)
        
      } else if (xi.method == "neural_network") {
        
        xi.models <- list()
        
        # For the lower bound
        y.train.down1 <- scale(train.valid.data1$y.tensor.train * gamma**(-sign(train.valid.data1$y.tensor.train - nn.init$Q.predict1[train.valid.data1$train.ind, 1])))
        y.valid.down1 <- scale(train.valid.data1$y.tensor.valid * gamma**(-sign(train.valid.data1$y.tensor.valid - nn.init$Q.predict1[train.valid.data1$valid.ind, 1])))
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.down1 <- train.nn(nn.architecture=base_neural_network,
                                           X.tensor.train=train.valid.data1$X.tensor.train,
                                           X.tensor.valid=train.valid.data1$X.tensor.valid,
                                           t.tensor.train=train.valid.data1$t.tensor.train,
                                           t.tensor.valid=train.valid.data1$t.tensor.valid,
                                           y.tensor.train=torch_tensor(y.train.down1, device=device),
                                           y.tensor.valid=torch_tensor(y.valid.down1, device=device),
                                           max.iter=nn.params$max.iter.resp,
                                           patience=nn.params$patience.resp,
                                           K=optimal.params.resp$K.optim,
                                           lr=optimal.params.resp$lr.optim,
                                           dim.hidden=optimal.params.resp$dim.hidden.optim,
                                           device=device,
                                           verbose=verbose)
        
        xi.models$xi.gmm.down1 <- trained.xi.model.down1$gmm
        xi.models$y.train.down1 <- y.train.down1
        xi.models$y.valid.down1 <- y.valid.down1
        
        y.train.down2 <- scale(train.valid.data2$y.tensor.train * gamma**(-sign(train.valid.data2$y.tensor.train - nn.init$Q.predict2[train.valid.data2$train.ind, 1])))
        y.valid.down2 <- scale(train.valid.data2$y.tensor.valid * gamma**(-sign(train.valid.data2$y.tensor.valid - nn.init$Q.predict2[train.valid.data2$valid.ind, 1])))
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.down2 <- train.nn(nn.architecture=base_neural_network,
                                           X.tensor.train=train.valid.data2$X.tensor.train,
                                           X.tensor.valid=train.valid.data2$X.tensor.valid,
                                           t.tensor.train=train.valid.data2$t.tensor.train,
                                           t.tensor.valid=train.valid.data2$t.tensor.valid,
                                           y.tensor.train=torch_tensor(y.train.down2, device=device),
                                           y.tensor.valid=torch_tensor(y.valid.down2, device=device),
                                           max.iter=nn.params$max.iter.resp,
                                           patience=nn.params$patience.resp,
                                           K=optimal.params.resp$K.optim,
                                           lr=optimal.params.resp$lr.optim,
                                           dim.hidden=optimal.params.resp$dim.hidden.optim,
                                           device=device,
                                           verbose=verbose)
        
        xi.models$xi.gmm.down2 <- trained.xi.model.down2$gmm
        xi.models$y.train.down2 <- y.train.down2
        xi.models$y.valid.down2 <- y.valid.down2
        
        # For the upper bound
        y.train.up1 <- scale(train.valid.data1$y.tensor.train * gamma**(sign(train.valid.data1$y.tensor.train - nn.init$Q.predict1[train.valid.data1$train.ind, 2])))
        y.valid.up1 <- scale(train.valid.data1$y.tensor.valid * gamma**(sign(train.valid.data1$y.tensor.valid - nn.init$Q.predict1[train.valid.data1$valid.ind, 2])))
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.up1 <- train.nn(nn.architecture=base_neural_network,
                                         X.tensor.train=train.valid.data1$X.tensor.train,
                                         X.tensor.valid=train.valid.data1$X.tensor.valid,
                                         t.tensor.train=train.valid.data1$t.tensor.train,
                                         t.tensor.valid=train.valid.data1$t.tensor.valid,
                                         y.tensor.train=torch_tensor(y.train.up1, device=device),
                                         y.tensor.valid=torch_tensor(y.valid.up1, device=device),
                                         max.iter=nn.params$max.iter.resp,
                                         patience=nn.params$patience.resp,
                                         K=optimal.params.resp$K.optim,
                                         lr=optimal.params.resp$lr.optim,
                                         dim.hidden=optimal.params.resp$dim.hidden.optim,
                                         device=device,
                                         verbose=verbose)
        
        xi.models$xi.gmm.up1 <- trained.xi.model.up1$gmm
        xi.models$y.train.up1 <- y.train.up1
        xi.models$y.valid.up1 <- y.valid.up1
        
        y.train.up2 <- scale(train.valid.data2$y.tensor.train * gamma**(sign(train.valid.data2$y.tensor.train - nn.init$Q.predict2[train.valid.data2$train.ind, 2])))
        y.valid.up2 <- scale(train.valid.data2$y.tensor.valid * gamma**(sign(train.valid.data2$y.tensor.valid - nn.init$Q.predict2[train.valid.data2$valid.ind, 2])))
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.up2 <- train.nn(nn.architecture=base_neural_network,
                                           X.tensor.train=train.valid.data2$X.tensor.train,
                                           X.tensor.valid=train.valid.data2$X.tensor.valid,
                                           t.tensor.train=train.valid.data2$t.tensor.train,
                                           t.tensor.valid=train.valid.data2$t.tensor.valid,
                                           y.tensor.train=torch_tensor(y.train.up2, device=device),
                                           y.tensor.valid=torch_tensor(y.valid.up2, device=device),
                                           max.iter=nn.params$max.iter.resp,
                                           patience=nn.params$patience.resp,
                                           K=optimal.params.resp$K.optim,
                                           lr=optimal.params.resp$lr.optim,
                                           dim.hidden=optimal.params.resp$dim.hidden.optim,
                                           device=device,
                                           verbose=verbose)
        
        xi.models$xi.gmm.up2 <- trained.xi.model.up2$gmm
        xi.models$y.train.up2 <- y.train.up2
        xi.models$y.valid.up2 <- y.valid.up2
        
      } else {
        stop('xi.method must be "regression_forest" or "neural_network"')
      }
      
    }
    
    end.time.qb2 <- Sys.time()
    
  } else {  # For the bootstrap resample
    
    start.time.qb2 <- Sys.time()
    start.time.jesson2 <- Sys.time()
    
    # Put the networks in train mode
    nn.init$gps.gmm1$train()
    nn.init$gps.gmm2$train()
    nn.init$resp.gmm1$train()
    nn.init$resp.gmm2$train()
    
    optimal.params.gps <- NULL
    optimal.params.resp <- NULL
    
    # Retrain the response model
    # On D1
    trained.resp.model1 <- train.nn(nn.architecture=base_neural_network,
                                    nn.model=nn.init$resp.gmm1,
                                    X.tensor.train=train.valid.data1$X.tensor.train,
                                    X.tensor.valid=train.valid.data1$X.tensor.valid,
                                    t.tensor.train=train.valid.data1$t.tensor.train,
                                    t.tensor.valid=train.valid.data1$t.tensor.valid,
                                    y.tensor.train=train.valid.data1$y.tensor.train,
                                    y.tensor.valid=train.valid.data1$y.tensor.valid,
                                    max.iter=nn.params$max.iter.resp,
                                    patience=nn.params$patience.resp,
                                    K=nn.params$K.resp,
                                    lr=nn.params$lr.resp,
                                    dim.hidden=nn.params$hid.dim.resp,
                                    device=device,
                                    verbose=verbose)
    
    nn.init$resp.gmm1 <- trained.resp.model1$gmm
    
    # On D2
    trained.resp.model2 <- train.nn(nn.architecture=base_neural_network,
                                    nn.model=nn.init$resp.gmm2,
                                    X.tensor.train=train.valid.data2$X.tensor.train,
                                    X.tensor.valid=train.valid.data2$X.tensor.valid,
                                    t.tensor.train=train.valid.data2$t.tensor.train,
                                    t.tensor.valid=train.valid.data2$t.tensor.valid,
                                    y.tensor.train=train.valid.data2$y.tensor.train,
                                    y.tensor.valid=train.valid.data2$y.tensor.valid,
                                    max.iter=nn.params$max.iter.resp,
                                    patience=nn.params$patience.resp,
                                    K=nn.params$K.resp,
                                    lr=nn.params$lr.resp,
                                    dim.hidden=nn.params$hid.dim.resp,
                                    device=device,
                                    verbose=verbose)
    
    nn.init$resp.gmm2 <- trained.resp.model2$gmm
    
    end.time.jesson2 <- Sys.time()
    
    # Retrain the GPS model
    # On D1
    trained.gps.model1 <- train.nn(nn.architecture=base_neural_network_gps,
                                   nn.model=nn.init$gps.gmm1,
                                   X.tensor.train=train.valid.data1$X.tensor.train,
                                   X.tensor.valid=train.valid.data1$X.tensor.valid,
                                   t.tensor.train=train.valid.data1$t.tensor.train,
                                   t.tensor.valid=train.valid.data1$t.tensor.valid,
                                   max.iter=nn.params$max.iter.gps,
                                   patience=nn.params$patience.gps,
                                   K=nn.params$K.gps,
                                   lr=nn.params$lr.gps,
                                   dim.hidden=nn.params$hid.dim.gps,
                                   device=device,
                                   verbose=verbose)
    
    nn.init$gps.gmm1 <- trained.gps.model1$gmm
    
    # On D2
    trained.gps.model2 <- train.nn(nn.architecture=base_neural_network_gps,
                                   nn.model=nn.init$gps.gmm2,
                                   X.tensor.train=train.valid.data2$X.tensor.train,
                                   X.tensor.valid=train.valid.data2$X.tensor.valid,
                                   t.tensor.train=train.valid.data2$t.tensor.train,
                                   t.tensor.valid=train.valid.data2$t.tensor.valid,
                                   max.iter=nn.params$max.iter.gps,
                                   patience=nn.params$patience.gps,
                                   K=nn.params$K.gps,
                                   lr=nn.params$lr.gps,
                                   dim.hidden=nn.params$hid.dim.gps,
                                   device=device,
                                   verbose=verbose)
    
    nn.init$gps.gmm2 <- trained.gps.model2$gmm
    
    if (compute.QB) {

      # Split the conditional quantiles between D1 and D2
      nn.init$Q.predict1 <- Q.predict[data1.ind, ]
      nn.init$Q.predict2 <- Q.predict[data2.ind, ]
      
      
      if (is.null(xi.method)) {
        
        xi.models <- NULL
        
      } else if (xi.method == "regression_forest") {
        
        # Train function xi that brings double-robustness
        # Model fitting on D1
        X.t.d1 <- cbind(X1, t1)
        names(X.t.d1)[names(X.t.d1) == "t1"] <- "t"
        X.t.d1 <- as.matrix(X.t.d1)
        # For the lower bound
        Y.gamma.down.d1 <- scale(Y1 * gamma**(-sign(Y1 - nn.init$Q.predict1[, 1])))
        xi.model.down.d1 <- regression_forest(X=X.t.d1,
                                              Y=Y.gamma.down.d1,
                                              num.trees=100,
                                              tune.parameters="all")
        # For the upper bound
        Y.gamma.up.d1 <- scale(Y1 * gamma**(sign(Y1 - nn.init$Q.predict1[, 2])))
        xi.model.up.d1 <- regression_forest(X=X.t.d1,
                                            Y=Y.gamma.up.d1,
                                            num.trees=100,
                                            tune.parameters="all")
        
        # Model fitting on D2
        X.t.d2 <- cbind(X2, t2)
        names(X.t.d2)[names(X.t.d2) == "t2"] <- "t"
        X.t.d2 <- as.matrix(X.t.d2)
        # For the lower bound
        Y.gamma.down.d2 <- scale(Y2 * gamma**(-sign(Y2 - nn.init$Q.predict2[, 1])))
        xi.model.down.d2 <- regression_forest(X=X.t.d2,
                                              Y=Y.gamma.down.d2,
                                              num.trees=100,
                                              tune.parameters="all")
        # For the upper bound
        Y.gamma.up.d2 <- scale(Y2 * gamma**(sign(Y2 - nn.init$Q.predict2[, 2])))
        xi.model.up.d2 <- regression_forest(X=X.t.d2,
                                            Y=Y.gamma.up.d2,
                                            num.trees=100,
                                            tune.parameters="all")
        
        xi.models <- list(xi.model.down.d1=xi.model.down.d1,
                          xi.model.up.d1=xi.model.up.d1,
                          xi.model.down.d2=xi.model.down.d2,
                          xi.model.up.d2=xi.model.up.d2,
                          Y.gamma.down.d1=Y.gamma.down.d1,
                          Y.gamma.up.d1=Y.gamma.up.d1,
                          Y.gamma.down.d2=Y.gamma.down.d2,
                          Y.gamma.up.d2=Y.gamma.up.d2)
        
      } else if (xi.method == "neural_network") {
        
        xi.models$xi.gmm.down1$train()
        xi.models$xi.gmm.down2$train()
        xi.models$xi.gmm.up1$train()
        xi.models$xi.gmm.up2$train()
        
        # For the lower bound
        y.train.down1 <- scale(train.valid.data1$y.tensor.train * gamma**(-sign(train.valid.data1$y.tensor.train - nn.init$Q.predict1[train.valid.data1$train.ind, 1])))
        y.valid.down1 <- scale(train.valid.data1$y.tensor.valid * gamma**(-sign(train.valid.data1$y.tensor.valid - nn.init$Q.predict1[train.valid.data1$valid.ind, 1])))
        y.train.down2 <- scale(train.valid.data2$y.tensor.train * gamma**(-sign(train.valid.data2$y.tensor.train - nn.init$Q.predict2[train.valid.data2$train.ind, 1])))
        y.valid.down2 <- scale(train.valid.data2$y.tensor.valid * gamma**(-sign(train.valid.data2$y.tensor.valid - nn.init$Q.predict2[train.valid.data2$valid.ind, 1])))
        # For the upper bound
        y.train.up1 <- scale(train.valid.data1$y.tensor.train * gamma**(sign(train.valid.data1$y.tensor.train - nn.init$Q.predict1[train.valid.data1$train.ind, 2])))
        y.valid.up1 <- scale(train.valid.data1$y.tensor.valid * gamma**(sign(train.valid.data1$y.tensor.valid - nn.init$Q.predict1[train.valid.data1$valid.ind, 2])))
        y.train.up2 <- scale(train.valid.data2$y.tensor.train * gamma**(sign(train.valid.data2$y.tensor.train - nn.init$Q.predict2[train.valid.data2$train.ind, 2])))
        y.valid.up2 <- scale(train.valid.data2$y.tensor.valid * gamma**(sign(train.valid.data2$y.tensor.valid - nn.init$Q.predict2[train.valid.data2$valid.ind, 2])))
        
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.down1 <- train.nn(nn.architecture=base_neural_network,
                                           nn.model=xi.models$xi.gmm.down1,
                                           X.tensor.train=train.valid.data1$X.tensor.train,
                                           X.tensor.valid=train.valid.data1$X.tensor.valid,
                                           t.tensor.train=train.valid.data1$t.tensor.train,
                                           t.tensor.valid=train.valid.data1$t.tensor.valid,
                                           y.tensor.train=torch_tensor(y.train.down1, device=device),
                                           y.tensor.valid=torch_tensor(y.valid.down1, device=device),
                                           max.iter=nn.params$max.iter.resp,
                                           patience=nn.params$patience.resp,
                                           K=nn.params$K.resp,
                                           lr=nn.params$lr.resp,
                                           dim.hidden=nn.params$hid.dim.resp,
                                           device=device,
                                           verbose=verbose)
        
        xi.models$xi.gmm.down1 <- trained.xi.model.down1$gmm
        xi.models$y.train.down1 <- y.train.down1
        xi.models$y.valid.down1 <- y.valid.down1
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.down2 <- train.nn(nn.architecture=base_neural_network,
                                           nn.model=xi.models$xi.gmm.down2,
                                           X.tensor.train=train.valid.data2$X.tensor.train,
                                           X.tensor.valid=train.valid.data2$X.tensor.valid,
                                           t.tensor.train=train.valid.data2$t.tensor.train,
                                           t.tensor.valid=train.valid.data2$t.tensor.valid,
                                           y.tensor.train=torch_tensor(y.train.down2, device=device),
                                           y.tensor.valid=torch_tensor(y.valid.down2, device=device),
                                           max.iter=nn.params$max.iter.resp,
                                           patience=nn.params$patience.resp,
                                           K=nn.params$K.resp,
                                           lr=nn.params$lr.resp,
                                           dim.hidden=nn.params$hid.dim.resp,
                                           device=device,
                                           verbose=verbose)
        
        xi.models$xi.gmm.down2 <- trained.xi.model.down2$gmm
        xi.models$y.train.down2 <- y.train.down2
        xi.models$y.valid.down2 <- y.valid.down2
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.up1 <- train.nn(nn.architecture=base_neural_network,
                                         nn.model=xi.models$xi.gmm.up1,
                                         X.tensor.train=train.valid.data1$X.tensor.train,
                                         X.tensor.valid=train.valid.data1$X.tensor.valid,
                                         t.tensor.train=train.valid.data1$t.tensor.train,
                                         t.tensor.valid=train.valid.data1$t.tensor.valid,
                                         y.tensor.train=torch_tensor(y.train.up1, device=device),
                                         y.tensor.valid=torch_tensor(y.valid.up1, device=device),
                                         max.iter=nn.params$max.iter.resp,
                                         patience=nn.params$patience.resp,
                                         K=nn.params$K.resp,
                                         lr=nn.params$lr.resp,
                                         dim.hidden=nn.params$hid.dim.resp,
                                         device=device,
                                         verbose=verbose)
        
        xi.models$xi.gmm.up1 <- trained.xi.model.up1$gmm
        xi.models$y.train.up1 <- y.train.up1
        xi.models$y.valid.up1 <- y.valid.up1
        
        # Train the final model on 90% of D1 and validate on 10% of D1
        trained.xi.model.up2 <- train.nn(nn.architecture=base_neural_network,
                                         nn.model=xi.models$xi.gmm.up2,
                                         X.tensor.train=train.valid.data2$X.tensor.train,
                                         X.tensor.valid=train.valid.data2$X.tensor.valid,
                                         t.tensor.train=train.valid.data2$t.tensor.train,
                                         t.tensor.valid=train.valid.data2$t.tensor.valid,
                                         y.tensor.train=torch_tensor(y.train.up2, device=device),
                                         y.tensor.valid=torch_tensor(y.valid.up2, device=device),
                                         max.iter=nn.params$max.iter.resp,
                                         patience=nn.params$patience.resp,
                                         K=nn.params$K.resp,
                                         lr=nn.params$lr.resp,
                                         dim.hidden=nn.params$hid.dim.resp,
                                         device=device,
                                         verbose=verbose)
        
        xi.models$xi.gmm.up2 <- trained.xi.model.up2$gmm
        xi.models$y.train.up2 <- y.train.up2
        xi.models$y.valid.up2 <- y.valid.up2
        
      } else {
        stop('xi.method must be "regression_forest" or "neural_network"')
      }

    }
    
    end.time.qb2 <- Sys.time()
    
  }
  
  # Put the networks in evaluation mode
  nn.init$gps.gmm1$eval()
  nn.init$gps.gmm2$eval()
  nn.init$resp.gmm1$eval()
  nn.init$resp.gmm2$eval()
  
  if (!is.null(xi.method)) {
    if (xi.method == "neural_network") {
      xi.models$xi.gmm.down1$eval()
      xi.models$xi.gmm.down2$eval()
      xi.models$xi.gmm.up1$eval()
      xi.models$xi.gmm.up2$eval()
    }
  }
  
  start.time.qb3 <- Sys.time()
  
  if (compute.QB) {
    
    # To save the estimation of the APO under unconfoundedness
    qb.unconf.per.window.and.dose <- matrix(NA, nrow=length(bandwidths), ncol=2*doses.length)

    # Matrix of size n.windows * (2*doses.lengths)
    qb.PEI.per.window.and.dose <- foreach(d=1:doses.length, .combine="cbind") %do% {
      
      # Compute the PEI for each given bandwidth
      qb.PEI.per.window <- APO.QB.PEI.per.window(X=X,
                                                 Y=Y,
                                                 t=t,
                                                 dose=doses[d],
                                                 alpha=NULL,  # Useless here
                                                 lambda=gamma,
                                                 window=NULL,  # Useless here
                                                 windows.range=bandwidths,
                                                 cond.quant.method=cond.quant.method,
                                                 stabilization=stabilization,
                                                 use.AIPW=TRUE,
                                                 AIPW.outcome.algo="nn_GMM",
                                                 cond.dens.method="nn_GMM",
                                                 nn.init=nn.init,
                                                 xi.method=xi.method,
                                                 xi.models=xi.models,
                                                 d1.d2.ind=d1.d2.ind,
                                                 K=NULL,
                                                 compute.CI=NULL,
                                                 B=NULL,
                                                 use.parallel=NULL,
                                                 e.variables=NULL,
                                                 use.cross.fitting=NULL,
                                                 device=device,
                                                 verbose=verbose)
      
      qb.unconf.per.window.and.dose[, (2*d-1):(2*d)] <- qb.PEI.per.window$PEI.bounds.per.window[, 3:4]
      
      # Matrix of size n.windows*2
      qb.PEI.per.window$PEI.bounds.per.window[, 1:2]
    }
    
    qb.PEI.per.window.and.dose.lb <- qb.PEI.per.window.and.dose[, seq(from=1, to=ncol(qb.PEI.per.window.and.dose), by=2)]
    qb.PEI.per.window.and.dose.ub <- qb.PEI.per.window.and.dose[, seq(from=2, to=ncol(qb.PEI.per.window.and.dose), by=2)]
    qb.unconf.per.window.and.dose.lb <- qb.unconf.per.window.and.dose[, seq(from=1, to=ncol(qb.unconf.per.window.and.dose), by=2)]
    qb.unconf.per.window.and.dose.ub <- qb.unconf.per.window.and.dose[, seq(from=2, to=ncol(qb.unconf.per.window.and.dose), by=2)]
    
  } else {
    qb.PEI.per.window.and.dose.lb <- matrix(0, nrow=length(bandwidths), ncol=doses.length)
    qb.PEI.per.window.and.dose.ub <- matrix(0, nrow=length(bandwidths), ncol=doses.length)
    qb.unconf.per.window.and.dose.lb <- matrix(0, nrow=length(bandwidths), ncol=doses.length)
    qb.unconf.per.window.and.dose.ub <- matrix(0, nrow=length(bandwidths), ncol=doses.length)
  }
  
  end.time.qb3 <- Sys.time()
  
  start.time.jesson3 <- Sys.time()
  
  # If we want to compute bounds for Jesson et al.'s method
  if (compute.Jesson) {
    
    jesson.unconf.per.dose <- rep(NA, doses.length)
    
    # Vector of size 2*doses.length
    jesson.PEI.per.dose <- foreach(d=1:doses.length, .combine="c") %do% {
      
      jesson.PEI <- APO.Jesson.PEI(dose=doses[d],
                                   lambda=gamma,
                                   capo.method="nn_gmm",
                                   nn.init=nn.init,
                                   d1.d2.ind=d1.d2.ind,
                                   unconf.capo.model=NULL,
                                   X.tensor=X.tensor,
                                   X.sample=NULL,  # X.sample is only useful for regression_forest
                                   Y.sample.len=Y.sample.len,
                                   use.parallel=use.parallel.jesson,
                                   device=device,
                                   verbose=verbose)
      
      jesson.unconf.per.dose[d] <- jesson.PEI$unconf.APO.estimate
      
      jesson.PEI$PEI.bounds
    }
    
  } else {
    jesson.unconf.per.dose <- rep(0, doses.length)
    jesson.PEI.per.dose <- rep(0, 2*doses.length)
  }
  
  jesson.PEI.per.dose.lb <- jesson.PEI.per.dose[seq(from=1, to=length(jesson.PEI.per.dose), by=2)]
  jesson.PEI.per.dose.ub <- jesson.PEI.per.dose[seq(from=2, to=length(jesson.PEI.per.dose), by=2)]
  
  end.time.jesson3 <- Sys.time()
  
  exec.time.qb <- as.numeric(end.time.qb3 - start.time.qb3, units="secs") + as.numeric(end.time.qb2 - start.time.qb2, units="secs") + as.numeric(end.time.qb1 - start.time.qb1, units="secs")
  exec.time.jesson <- as.numeric(end.time.jesson3 - start.time.jesson3, units="secs") + as.numeric(end.time.jesson2 - start.time.jesson2, units="secs") + as.numeric(end.time.jesson1 - start.time.jesson1, units="secs")
  
  return(list(qb.PEI.per.window.and.dose.lb=qb.PEI.per.window.and.dose.lb,
              qb.PEI.per.window.and.dose.ub=qb.PEI.per.window.and.dose.ub,
              qb.unconf.per.window.and.dose.lb=qb.unconf.per.window.and.dose.lb,
              qb.unconf.per.window.and.dose.ub=qb.unconf.per.window.and.dose.ub,
              jesson.PEI.per.dose.lb=jesson.PEI.per.dose.lb,
              jesson.PEI.per.dose.ub=jesson.PEI.per.dose.ub,
              jesson.unconf.per.dose=jesson.unconf.per.dose,
              X2=X2, t2=t2, Y2=Y2,
              d1.d2.ind=d1.d2.ind,
              Q.predict=Q.predict,
              nn.init=nn.init,
              xi.models=xi.models,
              optimal.params.gps=optimal.params.gps,
              optimal.params.resp=optimal.params.resp,
              exec.time.qb=exec.time.qb,
              exec.time.jesson=exec.time.jesson))
}


train.valid.split <- function(n.data, X.tensor, t.tensor, y.tensor) {

  # Divide D1 into train and validation sets
  n.train <- floor(0.9*n.data)
  n.valid <- n.data - n.train

  train.ind <- sort(sample(1:n.data, n.train))
  valid.ind <- setdiff(1:n.data, train.ind)

  X.tensor.train <- X.tensor[train.ind, ]
  X.tensor.valid <- X.tensor[valid.ind, ]

  t.tensor.train <- t.tensor[train.ind]
  t.tensor.valid <- t.tensor[valid.ind]

  y.tensor.train <- y.tensor[train.ind]
  y.tensor.valid <- y.tensor[valid.ind]

  return(list(X.tensor.train=X.tensor.train, X.tensor.valid=X.tensor.valid,
              t.tensor.train=t.tensor.train, t.tensor.valid=t.tensor.valid,
              y.tensor.train=y.tensor.train, y.tensor.valid=y.tensor.valid,
              train.ind=train.ind, valid.ind=valid.ind))
}

