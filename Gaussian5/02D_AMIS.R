# Fit model with importance distribution parameters learnt from data

# Load libraries
library("INLA")
library("mvtnorm")

# Load AMIS functions
source("../genFuncs.R")
source("../inlaAMIS.R")

# Load data
load("data.RData")


# Define prior for parameters
prior <- function(x, log = TRUE) {
  if(log == TRUE) {
    res <- 0
  } else {
    res <- 1
  }
  return(res)
}

# Define proposal distribution for log(tau_i)
# Density
dprop <- function(y, x, sigma = diag(5, 2, 2), log = TRUE) {
  return(mvtnorm::dmvnorm(y, mean = x, sigma = sigma, log = log))
}

# Sampling random values
rprop <- function(x, sigma = diag(5, 2, 2)){
  return(mvtnorm::rmvnorm(1, mean=x, sigma = sigma))
}



# Fot model with INLA given gamma
# g_values: Vector of values for gamma

# Data in matrix form per groups
Y <- matrix(NA, ncol = n_groups, nrow = n)
offset <- 0
for(i in 1:n_groups) {
  Y[offset + 1:(ni[i]), i] <- y[offset + 1:(ni[i])]
  offset <- offset + ni[i]
}

d <- as.list(d)
d$Y <- Y
d$x_prec1 <- x_prec1

# g_values: log(tau_1), ..., log(tau_ngroup)
fit_model <- function(data, g_values) {
  require("INLA")
  # Values of gamma_0 and gamma_1
  #g_values <- samples[1, ]

  # Prior for the precision 
  hyper_prec <- lapply(g_values, function(X) {
    list(hyper = list(prec = list(initial = X, fixed = TRUE)))
  })

  # Fit INLA model (y_ij)
  m1 <- inla(Y ~ 1 + x_mean1,
    data = list(Y = data$Y, x_mean1 = data$x_mean1, grp = data$grp),
    family = rep("gaussian", n_groups),
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    control.family = hyper_prec)
  #print("M1 fitted")

  # Fit model (log(tau_i))
  m2 <- inla(g_values ~ 1 + x_prec1 + f(idx, model = "iid"),
    data = list(g_values = as.vector(g_values), x_prec1 = data$x_prec1, idx = 1:n_groups),
    family = "gaussian",
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    control.family = list(hyper = list(prec = list(initial = 15, fixed = TRUE)))
  )
  #print("M2 fitted")

  res <- list(
    mlik = m1$mlik[1, 1] + m2$mlik[1, 1],
    dists = c(m1$marginals.fixed, m2$marginals.fixed, m2$marginals.hyperpar,
      m2$internal.marginals.hyperpar)
  )
}

# Precisions versus estimated precisions
g_precs_hat <- 1 / as.vector(with(d, by(y, grp, var)))
plot(g_precs, g_precs_hat)
abline(0, 1)

# Fit model using AMIS
init <- list(
  # Vague importance distribution
  #mu = rep(0, n_groups),
  #cov = diag(5, n_groups, n_groups)
  # Use estimated log(precs)
  mu = log(g_precs_hat),
  cov = diag(var(log(g_precs_hat)) / ni[1], n_groups, n_groups)
)

# Number of samples
N0 <- 1000 #Initial step
Nt <- rep(1000, 10) #Adaptations

set.seed(123)
res_time <- system.time(
  res <- inlaAMIS(data = d, init = init, prior = prior, d.prop = dprop,
    r.prop = rprop, fit.inla = fit_model, N_t = Nt,
    N_0 = N0, ncores = 60)
)

save(file = "AMISD.RData", list = ls())

