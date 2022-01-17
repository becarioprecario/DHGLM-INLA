# Summary of results for paper

# Load results
load("AMIS_PO.RData")
load("JAGS_PO.RData")


# AMIS: Summary stats

# Computing times
res_time

# Number of second per simulation
res_time[3] / sum(N0, Nt)

#Re-scaled weights
ww <- res$weight / sum(res$weight)
summary(ww)

# Effective smaple size
sum(ww)^2  / sum(ww^2)

# Determinant of variance matrix of importance distribution
#  This may be useful to monitor convergence
pdf(file = "det_PO.pdf")
plot(apply(res$theta$a.cov, 3, det))
dev.off()

# Summary stats (mean and variance matrix)
res$theta$a.mu[nrow(res$theta$a.mu), ]
res$theta$a.cov[, , nrow(res$theta$a.mu)]

# Summary stats using inla.zmarginal
library("INLA")
# beta
lapply(res$margs, inla.zmarginal)
# gamma
lapply(res$eta_kern, inla.zmarginal)

# JAGS: Summary stats
# Post. means
js$beta
js$gamma

# Post. st. dev.
# beta
apply(js$beta[, , 1], 1, sd)
#gamma
apply(js$gamma[, , 1], 1, sd)

# Post. variance
var(t(js$gamma[, , 1]))

pdf(file = "marginals_PO.pdf")
# PLot the posterior densities
par(mfrow = c(2, 2))
# beta_0
plot(density(js$beta[1, , 1]))
lines(res$margs[[1]], col = "red")
#beta_1
plot(density(js$beta[2, , 1]))
lines(res$margs[[2]], col = "red")
# gamma_0
plot(density(js$gamma[1, , 1]))
lines(res$eta_kern[[1]], col = "red")
# gamma_1
plot(density(js$gamma[2,,1]))
lines(res$eta_kern[[2]], col = "red")

dev.off()



