# Load data
load("data.RData")

# JAGS
library(rjags)
library(coda)

# Define model
jm <- jags.model("../bugs_models/model_gaussian_reff.bug",
  data = list(y = y, x_mean1 = x_mean1, x_prec1 = x_prec1, n = n,
    ngroups = n_groups, grp = rep(1:n_groups, ni)),
  inits = list(beta = c(0, 0), gamma = c(0, 0),
    ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1)
)

# Burn-in
update(jm, n.iter = 10000)

# Sample from model
js <- jags.samples(jm,
  variable.names = c("beta", "gamma", "prec", "precu", "reff"),
  n.iter = 100000, thin = 100)

save(file = "JAGS.RData", list = ls())

