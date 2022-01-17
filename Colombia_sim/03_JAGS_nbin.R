# Load data
load("data.RData")

# JAGS
library(rjags)
library(coda)

# Define model
jm <- jags.model("../bugs_models/model_nbinom_offset.bug",
  data = list(y = y, x_mean1 = x_mean1, x_prec1 = x_prec1, n = n,
    log_offset = log(offset1)),
  inits = list(beta = c(0, 0), gamma = c(0, 0),
    ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1)
)

# Burn-in
update(jm, n.iter = 10000)

# Sample from model
js <- jags.samples(jm, variable.names = c("beta", "gamma", "size"),
  n.iter = 100000, thin = 100)

save(file = "JAGS_NB.RData", list = ls())

