# Fit model using AMIS with INLA

# Load libraries
library("INLA")
library("mvtnorm")

# Load AMIS functions
source("../genFuncs.R")
source("../inlaAMIS.R")

# Load data
load("data.RData")


#Subset data
#n_groups <- 8
#d <- subset(d, grp %in% levels(d$grp)[1:n_groups])
#d$grp <- as.factor(as.character(d$grp))

# Fit model using lme4
library("lme4")
mm4 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)

# Fit linear models to each subject (these estimates may be used
#  later)
lm_subject <- lapply(unique(d$grp), function(X) {
  aux <- subset(d, grp == X)
  m0 <- lm(y ~ x_mean1, data = aux)
  res <- c(coefficients(m0), -log(var(m0$residuals)), -log(var(aux$y)))
  return(res)
})

lm_subject <- do.call(rbind, lm_subject)




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


# Data in matrix form per groups
Y <- matrix(NA, ncol = n_groups, nrow = n)
offset <- 0
for(i in 1:n_groups) {
  Y[offset + 1:(ni[i]), i] <- y[offset + 1:(ni[i])]
  offset <- offset + ni[i]
}

d <- as.list(d)
d$Y <- Y

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
  m1 <- inla(Y ~ 1 + f(grp, x_mean1, model = "iid"),
    data = list(Y = data$Y, x_mean1 = data$x_mean1, grp = data$grp),
    family = rep("gaussian", n_groups),
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    control.family = hyper_prec)
  #print("M1 fitted")
  
  # IMPORTANT: INLA does not use this term when computing the marg. lik.
  #  because it is a constant term.
  # So it is added manually.
  #m1$logdet <- as.numeric(Matrix::determinant(Cmatrix)$modulus)
  #m1$mlik <- m1$mlik + m1$logdet / 2
  
  # Fit model (log(tau_i))
  m2 <- inla(g_values ~ 1 + f(idx, model = "iid",
      hyper = list(prec = list(param = c(0.01, 0.01)))),
    data = list(g_values = as.vector(g_values), idx = 1:n_groups),
    family = "gaussian",
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    control.family = list(hyper = list(prec = list(initial = 15, fixed = TRUE)))
  )
  #print("M2 fitted")
  
  res <- list(
    mlik = m1$mlik[1, 1] + m2$mlik[1, 1],
    dists = c(m1$marginals.fixed, m1$marginals.hyperpar, 
      m2$marginals.fixed, m2$marginals.hyperpar,
      m2$internal.marginals.hyperpar)
  )
}

#log_prec_hat <- 
#  as.vector(log(1 / by(sleepstudy$Reaction, sleepstudy$Subject, var)))

# Variance of residuals of lm
#log_prec_hat <- lm_subject[, 3]
# Variance of observed data
log_prec_hat <- lm_subject[, 4]


# Test some initial values for the mean by resampling log_prec_hat
set.seed(123)
init_lp <- c(list(log_prec_hat),
  lapply(1:500, function(X) {sample(log_prec_hat)})
)

# Compute mlik
init_lp_mlik <- mclapply(1:length(init_lp), function(X) {
  fit_model(d, init_lp[[X]])$mlik
  }, mc.cores = 60
)


# Params
init_mean <- init_lp[[which.max(unlist(init_lp_mlik))]]
# sd equal to 10% of actual value to allow for small jumps and tuning
init_vcov <- diag( (0.05 * init_mean)^2, n_groups, n_groups)

print(max(unlist(init_lp_mlik)))
print(init_mean)


# Fit model using AMIS
init <- list(
  #mu = rep(0, n_groups),
  #cov = diag(3, n_groups, n_groups)
  # Use data to propose parameters of sampling distribution
  #mu = log_prec_hat,
  #cov = diag(var(log_prec_hat) / ni, n_groups, n_groups)
  mu = init_mean,
  cov = init_vcov
)

# Number of samples
N0 <- 1000 #Initial step
Nt <- rep(1000, 20) #Adaptations

# # Number of samples I have tried
# N0 <- 50 #Initial step
# Nt <- rep(100, 3) #Adaptations

# Fit model
set.seed(1)
res_time <- system.time(
  res <- inlaAMIS(data = d, init = init, prior = prior, d.prop = dprop,
    r.prop = rprop, fit.inla = fit_model, N_t = Nt,
    N_0 = N0, ncores = 60)
)


save(file = "AMIS_rcoef.RData", list = ls())
