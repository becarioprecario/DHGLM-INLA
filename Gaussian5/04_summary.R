#
# Gaussian model with hierarchy on the precisions
#


# Load results
load("JAGS.RData")
amis <- lapply(paste("AMIS", c("", "B", "C", "D", "B10","B55K"), ".RData", sep = ""),
  function(X) {mget(load(X))})

# Transform marginals 'prec internal scale' -> 'prec.'
library("INLA")
for(i in 1:length(amis)) {
   amis[[i]]$res$margs[[5]]<- data.frame(inla.tmarginal(exp, amis[[i]]$res$margs[[6]]))

}

# AMIS: Summary stats

# Computing times
lapply(amis, function(X) {X$res_time})

# Number of second per simulation
lapply(amis, function(X) {
  with(X, res_time[3] / sum(N0, Nt))
})

#Re-scaled weights
lapply(amis, function(X) {
  with(X, summary(res$weight / sum(res$weight)))#; summary(ww))
})

# Compute ESS
lapply(amis, function(X) {
  ww <- with(X, (res$weight / sum(res$weight)))#; summary(ww))
  print(sum(ww)^2 / sum(ww^2))
})


# Determinant of variance matrix of importance distribution
#  This may be useful to monitor convergence
# There are actual zeros for cases A and C
pdf(file = "det.pdf")
plot(log(apply(amis[[1]]$res$theta$a.cov, 3, det)), type = "n")
lapply(1:4, function(X) {
  lines(log(apply(amis[[X]]$res$theta$a.cov, 3, det)), lty = X)
 })
dev.off()


# Summary stats (mean and variance matrix)
lapply(amis, function(X) {
  with(X, res$theta$a.mu[nrow(res$theta$a.mu), ])
  with(X, res$theta$a.cov[, , nrow(res$theta$a.mu)])
})

# Summary stats using inla.zmarginal
library("INLA")

# beta, gamma and tau_u
lapply(amis, function(X) {
  with(X, lapply(res$margs, inla.zmarginal))
})
# log(prec_i)
lapply(amis, function(X) {
  with(X, lapply(res$eta_kern, inla.zmarginal))
})

# Compute marginals of prec_i
margs_prec <- lapply(amis, function(X) {
  with(X, lapply(res$eta_kern,
    function(X) {
      inla.tmarginal(exp, X)
    }
  ))
})

# Compute summary statistics of prec_imargss, inla.zmarginal)
lapply(margs_prec, function(X) {
  lapply(X, inla.zmarginal)
})


# JAGS: Summary stats
# Summary stats
js$beta
js$gamma
js$precu
js$reff
js$prec

# Variance-covariance
apply(js$gamma[,,1], 1, var)
cov(t(js$gamma[,,1]))

# Plot estimates of random effects
plot(g_reff, apply(js$reff[,,1], 1, mean))
abline(0,1)


# Post. st. dev.
# beta
apply(js$beta[, , 1], 1, sd)
#gamma
apply(js$gamma[, , 1], 1, sd)

# Post. variance
var(t(js$gamma[, , 1]))

pdf(file = "marginals.pdf")
# PLot the posterior densities of log(prec_i)
par(mfrow = c(2, 2))
for(i in 1:4) { #n_groups) {
  plot(density(js$prec[i, , 1]), type = "l")
  lapply(1:4, function(X) {
    lines((margs_prec[[X]])[[i]], col  ="red", lty = X)
  })
}
dev.off()


# Fixed effects + prec_u
pdf(file = "marginals_fixed.pdf")
par(mfrow = c(3, 2))
plot(density(js$beta[1, , 1]))
lapply(1:4, function(X) {
  lines(amis[[X]]$res$margs[[1]], col = "red", lty = X)
})
#beta_1
plot(density(js$beta[2, , 1]))
lapply(1:4, function(X) {
  lines(amis[[X]]$res$margs[[2]], col = "red", lty = X)
})
# gamma_0
plot(density(js$gamma[1, , 1]))
lapply(1:4, function(X) {
  lines(amis[[X]]$res$margs[[3]], col = "red", lty = X)
})
# gamma_1
plot(density(js$gamma[2,,1]))
lapply(1:4, function(X) {
  lines(amis[[X]]$res$margs[[4]], col = "red", lty = X)
})
# prec_u
plot(density(js$precu[ , ,1]))
lapply(1:4, function(X) {
  lines(amis[[X]]$res$margs[[5]], col = "red", lty = X)
})
dev.off()

