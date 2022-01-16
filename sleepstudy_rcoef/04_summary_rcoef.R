# Summary of results for paper

# Load results
load("AMIS_rcoef.RData")
load("JAGS_rcoef.RData")

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
pdf(file = "det_rcoef.pdf")
plot(apply(res$theta$a.cov, 3, det))
dev.off()

# Summary stats (mean and variance matrix)
res$theta$a.mu[nrow(res$theta$a.mu), ]
res$theta$a.cov[, , nrow(res$theta$a.mu)]

# Summary stats using inla.zmarginal
library("INLA")
# beta, precd, gamma, precu
lapply(res$margs, inla.zmarginal)
# delta (random slopes for subject)
lapply(res$eta_kern, inla.zmarginal)

# JAGS: Summary stats
# Post. means
js$beta
js$precd
js$gamma
js$precu


# MABEL: No he sido capaz de crear un objeto coda para reproducir este código
if(FALSE) {

# MCMC Summary stats with coda
library(coda)
#load("coda_sleepstudy_rcoef.Rda")
coda_sleepstudy_rcoef <- as.mcmc(do.call(cbind, js))

library(mcmcplots)
dev.new(height = 6, width = 9)
#windows(height = 6, width = 9)
traplot(coda_sleepstudy_rcoef, parms=c("beta", "gamma", "precu","precd"))
dev.new(height = 6, width = 9)
#windows(height = 6, width = 9)
denplot(coda_sleepstudy_rcoef, parms=c("beta", "gamma", "precu","precd"))

## Summary for fixed effects and precisions
summ.mcmc <- cbind(summary(coda_sleepstudy_rcoef)$statistics[c(1,21,2,22),c(1,2)],summary(coda_sleepstudy_rcoef)$quantiles[c(1,21,2,22),c(1,5)])
summ.mcmc

## Summary for the random coefficients for Days (delta)
summ.delta <- summary(coda_sleepstudy_rcoef)$statistics[23:40,c(1,2)]
summ.delta

} #if (FALSE)

# AMIS-INLA summary stats 
estimations <- lapply(res$margs, inla.zmarginal)
beta_est <- estimations[[1]]
precd_est <- estimations[[2]]
gamma_est <- estimations[[3]]
precu_est <- estimations[[4]]

reff_days <- lapply(res$eta_kern, inla.zmarginal)
delta_est <- reff_days 

## Summary for fixed effects and precisions
summ.is <- rbind(beta_est[c(1,2,3,7)],precd_est[c(1,2,3,7)],gamma_est[c(1,2,3,7)], precu_est[c(1,2,3,7)])
rownames(summ.is) <- c("beta","precd" , "gamma", "precu")
summ.is

## Summary for the random coefficients for Days (delta)
do.call("rbind", reff_days)[,c(1,2)]



# Post. st. dev.
# beta
apply(js$beta[, , 1, drop = FALSE], 1, sd)
# precd
apply(js$precd[, , 1, drop = FALSE], 1, sd)
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
plot(density(js$precd[, , 1]))
lines(res$margs[[2]], col = "red")
# gamma_0
plot(density(js$gamma[, , 1]))
lines(res$margs[[3]], col = "red")
# gamma_1
plot(density(js$precu[ , , 1]))
lines(res$margs[[4]], col = "red")

dev.off()


# Plot densities of log-precision
dev.new(height = 10, width = 10)
par(mfrow = c(5, 4))
for(i in 1:18) {
  plot(density(log(js$prec[i, , 1])), xaxt = 'n', yaxt = 'n', ann = FALSE)
  lines(res$eta_kern[[i]], col = "red")
}


