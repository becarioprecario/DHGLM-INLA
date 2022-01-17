# JAGS
library(rjags)
library(coda)

# Load data
load("data.RData")

# Define model
jm <- jags.model("../bugs_models/model_poisson.bug",
  data = list(y = y, x_mean1 = x_mean1, x_prec1 = x_prec1, n = n),
  inits = list(beta = c(0, 0), gamma = c(0, 0),
    ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1)
)

# Burn-in
update(jm, n.iter = 10000)

# Sample from model
js <- jags.samples(jm, variable.names = c("beta", "gamma", "prec"),
  n.iter = 100000, thin = 100)

# Summary stats
js$beta
js$gamma

apply(js$gamma[,,1], 1, var)

save(file = "JAGS.RData", list = ls())


