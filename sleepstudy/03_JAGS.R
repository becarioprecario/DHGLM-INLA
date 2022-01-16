# Fit model with JAGS

# Load libraries
library("rjags")

# Load data
load("data.RData")

n <- length(d$y)
ngroups <- nlevels(d$grp)

data_jags <- list(n = n, ngroups = ngroups, y = d$y,
  grp = d$grp, x_mean = d$x_mean1)


#parameters to observe
params <- c("beta", "gamma", "precu", "reff", "prec")

#initial values
inits <- list(beta = c(0, 0.01), gamma = 0, precu = 0.1,
  ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1)


jm <- jm <- jags.model("sleepstudy.bug",
  data = data_jags, inits = inits
)

# Burn-in
update(jm, n.iter = 10000)

# Samples for inference
js <- jags.samples(jm,
  variable.names = params,
  n.iter = 100000, n.thin = 100
)

js

save(file = "JAGS.RData", list = ls())
