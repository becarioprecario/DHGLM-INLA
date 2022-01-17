# Poisson model with random effects

# Load libraries
library(spdep)
library(INLA)

# Load AMIS functions
source("../genFuncs.R")
source("../inlaAMIS.R")

# Load data
load("data.RData")


# Define prior for parameters
prior <- function(x, log = TRUE) {
  if(log == TRUE) {
    res <- sum(dnorm(x, mean = 0, sd = sqrt(1 / 0.001), log = log))
  } else {
    res <- prod(dnorm(x, mean = 0, sd = sqrt(1 / 0.001), log = log))
  }
  return(res)
}


# Define proposal distribution
# Density
dprop <- function(y, x, sigma = diag(5, 2, 2), log = TRUE) {
  return(mvtnorm::dmvnorm(y, mean = x, sigma = sigma, log = log))
}

# Sampling random values
rprop <- function(x, sigma = diag(5,2,2)){
  return(mvtnorm::rmvnorm(1, mean=x, sigma = sigma))
}

# Fit INLA model
fit_model <- function(data, g_values) {
  # Values of gamma_0 and gamma_1
  #g_values <- samples[1, ]
  
  # Precision of random effect
  prec_reff <- exp(g_values[1] + g_values[2] * data$x_prec1)
  
  # Fit INLA model (no random effect)
  #m0 <- inla(y ~ 1 + x_mean1, data = d, family = "poisson")
  #summary(m0)
  
  
  # Fit INLA model (no random effect)
  # Precision matrix of random effects
  Cmatrix <- Diagonal(n, prec_reff)
  m1 <- inla(y ~ 1 + offset(log(offset1)) + x_mean1 +
      f(idx, model = "generic0", Cmatrix = Cmatrix,
    hyper = list(prec = list(initial = log(1), fixed = TRUE))),
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    data = data, family = "poisson")
  #summary(m1)
  
  # IMPORTANT: INLA does not use this term when computing the marg. lik.
  #  because it is a constant term.
  # So it is added manually.
  m1$logdet <- as.numeric(Matrix::determinant(Cmatrix)$modulus)
  m1$mlik <- m1$mlik + m1$logdet / 2
  
  res <- list(mlik = m1$mlik[1, 1],
              dists = m1$marginals.fixed
  )
  
  return(res)
}

# Fit model using AMIS
init <- list(mu = c(0, 0), 
  cov = diag(5, 2, 2)
)

# Number of samples
N0 <- 5000 #Initial step
Nt <- rep(1000, 10) #Adaptations

# Fit model
set.seed(1)
res_time <- system.time(
  res <- inlaAMIS(data = d, init = init, prior = prior, d.prop = dprop, 
    r.prop = rprop, fit.inla = fit_model, N_t = Nt,
    N_0 = N0, ncores = 60)
)

save(file = "AMIS_PO.RData", list = ls())
