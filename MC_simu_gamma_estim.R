# Import libraries
library(ggplot2)
library(latex2exp)  # For LaTeX expressions

# Import R files
source("./simulated_data_fun.R")


seed <- 1
set.seed(seed)

# Get provided arguments
args <- commandArgs()

# Whether to show or not messages and progress bars
verbose <- FALSE

# Number of Monte-Carlo samples
n.MC <- 1000

# Gamma used for each Monte-Carlo sample
gammas.per.parameter <- list()

parameter.to.test <- args[1] #"beta_U"  # "correlation.prop"

correlation.prop.vec <- c(0.05, 0.3, 0.5, 0.7, 0.95)
p_U <- 3
beta_U.list <- list(rep(0, p_U), rep(0.1, p_U), rep(0.2, p_U), rep(0.3, p_U))
#gamma_U.list <- list(rep(-0.75, p_U), rep(-0.5, p_U), rep(-0.25, p_U), rep(0, p_U), rep(0.25, p_U), rep(0.5, p_U), rep(0.75, p_U))

if (parameter.to.test == "beta_U") {
  parameter.to.change <- beta_U.list
} else if (parameter.to.test == "correlation.prop") {
  parameter.to.change <- correlation.prop.vec
} else {
  stop('parameter.to.test must be "beta_U" or "correlation.prop"')
}

n.param <- length(parameter.to.change)


for (i in 1:n.param) {
  
  message(paste("Parameter", i, "/", n.param))
  
  data.name <- "simul"
  
  if (verbose) {
    pb <- txtProgressBar(min=0,      # Minimum value of the progress bar
                         max=n.MC,   # Maximum value of the progress bar
                         style=3,    # Progress bar style (also available style = 1 and style = 2)
                         width=50,   # Progress bar width. Defaults to getOption("width")
                         char="=")   # Character used to create the bar
  }
  
  # Gamma used for each Monte-Carlo sample
  gammas <- rep(NA, n.MC)
  
  for (mc.ind in 1:n.MC) {
    
    ### Simulation
    
    # Sample size
    n <- 1000
    
    # Number of confounders
    p_X <- 5
    p_U <- 3
    
    # Correlation between confounders
    rho_X <- 0.3
    rho_U <- 0.3
    
    beta_X <- rep(0.3, p_X)
    
    if (parameter.to.test == "beta_U") {
      corr_XU_prop <- 0.5
      beta_U <- parameter.to.change[[i]]
    } else if (parameter.to.test == "correlation.prop") {
      corr_XU_prop <- parameter.to.change[i]
      beta_U <- rep(0.2, p_U)
    }
    
    gamma_X <- rep(0.2, p_X)
    gamma_U <- c(rep(0.4, floor(p_U/2)), rep(0.7, p_U-floor(p_U/2)))  # parameter.to.change[[i]]
    zeta <- -0.3
    # Observed Y
    sd_eps_T <- 0.5
    sd_eps_Y <- 0.3
    
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
    pdf_T_X <- dnorm(t, mean=mu_T_x(as.matrix(X), simu$beta_X, simu$beta_U, simu$cov_XU, simu$Sigma_X),
                     sd=sqrt(sigma_T_x(simu$sd_eps_T, simu$beta_U, simu$cov_XU, simu$Sigma_X, simu$Sigma_U))) + 1e-4
    
    quotient <- pdf_T_XU / pdf_T_X
    quant.order <- 0.99
    gamma.est <- unname(quantile(quotient, probs=c(quant.order)))
    gammas[mc.ind] <- gamma.est
    
    if (verbose) {
      setTxtProgressBar(pb, mc.ind)  # Add one unit to the progress bar
    }
    
  }
  
  gammas.per.parameter[[i]] <- gammas
  
  if (verbose) {
    close(pb)
  }
  
}

if (parameter.to.test == "beta_U") {
  
  gammas.df <- data.frame(gamma=gammas.per.parameter[[1]], param=parameter.to.change[[1]][1])
  
  for (i in 2:n.param) {
    gammas.df <- rbind(gammas.df, data.frame(gamma=gammas.per.parameter[[i]], param=parameter.to.change[[i]][1]))
  }
  
} else if (parameter.to.test == "correlation.prop") {
  
  gammas.df <- data.frame(gamma=gammas.per.parameter[[1]], param=parameter.to.change[1])
  
  for (i in 2:n.param) {
    gammas.df <- rbind(gammas.df, data.frame(gamma=gammas.per.parameter[[i]], param=parameter.to.change[i]))
  }
  
}


# Save data in a file
version <- parameter.to.test

gammas.per.param.file.name <- paste("./results/gammas_per_parameter", data.name, version, sep="_")
saveRDS(gammas.df, file=paste0(gammas.per.param.file.name, ".RData"))

our_method_blue <- "#00BFC4"

if (parameter.to.test == "beta_U") {
  
  # Plot with respect to beta_U
  labels.beta_u <- c()
  unique.param <- unique(gammas.df$param)
  for (i in 1:length(unique.param)) {
    labels.beta_u[i] <- paste("(", paste(rep(unique.param[i], p_U), collapse=", "), ")", sep="")
  }
  
  # png("./images/gamma_vs_beta_u.png", units="in", width=5.16, height=4.54, res=400)
  gamma.vs.beta_u.plot <- ggplot(gammas.df, aes(x=param, y=gamma, group=param)) +
    geom_boxplot(fill=our_method_blue) +
    scale_x_continuous(breaks=unique(gammas.df$param), labels=labels.beta_u) +
    scale_y_continuous(breaks=c(1, 5, 10, 15), labels=c(1, 5, 10, 15)) +
    labs(x=TeX(r"($\beta_U$)"), y=TeX(r"(Sensitivity parameter $\Gamma$)")) +
    theme_linedraw()
  gamma.vs.beta_u.plot
  # dev.off()
  
} else if (parameter.to.test == "correlation.prop") {
  
  # Plot with respect to the correlation value
  gammas.df$corr <- gammas.df$param * (1-rho_X)/p_U
  
  # png("./images/gamma_vs_correlation_xu.png", units="in", width=5.16, height=4.54, res=400)
  gamma.vs.corrXU.plot <- ggplot(gammas.df, aes(x=corr, y=gamma, group=corr)) +
    geom_boxplot(fill=our_method_blue) +
    scale_x_continuous(breaks=unique(gammas.df$corr),
                       labels=~paste(signif(unique(gammas.df$corr), 2), correlation.prop.vec, sep="\n"),
                       name=expression(atop(NA, atop(textstyle("Correlation"~rho[XU]), textstyle("Correlation parameter"~lambda))))) +  # Add another x label for the correlation proportion
    labs(x=TeX(r"(Correlation $\rho_{XU}$)"), y=TeX(r"(Sensitivity parameter $\Gamma$)")) +
    theme_linedraw() +
    theme(axis.title.x=element_text(vjust=4))
  gamma.vs.corrXU.plot
  # dev.off()
  
}


# # Plot with respect to gamma_U
# labels.gamma_u <- c()
# unique.param <- unique(gammas.df$param)
# for (i in 1:length(unique.param)) {
#   labels.gamma_u[i] <- paste("(", paste(rep(unique.param[i], p_U), collapse=", "), ")", sep="")
# }
# 
# png("images/gamma_vs_gamma_u.png", units="in", width=8.14, height=4.54, res=400)
# gamma.vs.gamma_u.plot <- ggplot(gammas.df, aes(x=param, y=gamma, group=param)) +
#   geom_boxplot(fill=our_method_blue) +
#   scale_x_continuous(breaks=unique(gammas.df$param), labels=labels.gamma_u) +
#   #scale_y_continuous(breaks=c(1, 5, 10, 15), labels=c(1, 5, 10, 15)) +
#   labs(x=TeX(r"($\gamma_U$)"), y=TeX(r"(Sensitivity parameter $\Gamma$)")) +
#   theme_linedraw()
# gamma.vs.gamma_u.plot
# dev.off()
