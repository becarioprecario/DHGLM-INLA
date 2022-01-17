# Fit model with JAGS

# Load libraries
library("rjags")

# Load data
load("data.RData")

n <- length(sleepstudy$Reaction)
ngroups <- nlevels(sleepstudy$Subject)

data_jags <- list(n = n, ngroups = ngroups, y = sleepstudy$Reaction,
  grp = sleepstudy$Subject, x_mean = sleepstudy$Days)


#parameters to observe
params <- c("beta", "delta", "gamma", "precu", "precd", "reff", "prec")

#initial values
inits <- list(beta = 0, gamma = 0, precu = 0.1, precd = 0.1,
  ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1)


jm <- jm <- jags.model("sleepstudy_rcoef.bug",
  data = data_jags, inits = inits
)

# Burn-in
update(jm, n.iter = 10000)

# Samples for inference
js <- jags.samples(jm,
  variable.names = params,
  n.iter = 100000, thin = 100
)

js

save(file = "JAGS_rcoef.RData", list = ls())
