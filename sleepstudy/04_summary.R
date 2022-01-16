# Summary of results for paper

# Load results
load("AMIS.RData")
load("JAGS.RData")

# Transform marginals 'prec internal scale' -> 'prec.'
library("INLA")
res$margs[[4]]<- data.frame(inla.tmarginal(exp, res$margs[[5]]))




# AMIS: Summary stats

# Computing times
res_time

# Number of second per simulation
res_time[3] / sum(N0, Nt)

#Re-scaled weights
ww <- res$weight / sum(res$weight)
summary(ww)

# Effective sample size 
sum(ww)^2 / sum(ww^2)


# Determinant of variance matrix of importance distribution
#  This may be useful to monitor convergence
pdf(file = "det.pdf")
plot(apply(res$theta$a.cov, 3, det))
dev.off()

# Summary stats (mean and variance matrix)
res$theta$a.mu[nrow(res$theta$a.mu), ]
res$theta$a.cov[, , nrow(res$theta$a.mu)]

# Summary stats using inla.zmarginal
library("INLA")
# beta, gamma, precu
lapply(res$margs, inla.zmarginal)
# log(prec)
lapply(res$eta_kern, inla.zmarginal)

# JAGS: Summary stats
# Post. means
js$beta
js$gamma
js$precu

# Post. st. dev.
# beta
apply(js$beta[, , 1], 1, sd)
#gamma
apply(js$gamma[, , 1, drop = FALSE], 1, sd)
# precu
apply(js$precu[, , 1, drop = FALSE], 1, sd)

pdf(file = "marginals.pdf")
# PLot the posterior densities
par(mfrow = c(2, 2))
# beta_0
plot(density(js$beta[1, , 1]))
lines(res$margs[[1]], col = "red")
#beta_1
plot(density(js$beta[2, , 1]))
lines(res$margs[[2]], col = "red")
# gamma_0
plot(density(js$gamma[, , 1]))
lines(res$margs[[3]], col = "red")
# gamma_1
plot(density(js$precu[ , , 1]))
lines(res$margs[[4]], col = "red")

dev.off()



