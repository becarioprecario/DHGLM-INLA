############################### Summary of results 
############################### Gaussian model with random slopes for each subject (Sleepstudy data)


library(INLA)
library(mcmcplots)
library(ggplot2)
library(coda)
library(rjags)
library(ggpubr) 


####### Load results
load("AMIS_rcoef.RData")
load("JAGS_rcoef.RData")


####### Computing times
res_time


####### Number of second per simulation
res_time[3] / sum(N0, Nt)


####### Re-scaled weights
ww <- res$weight / sum(res$weight)
summary(ww)


####### Determinant of variance matrix of importance distribution
#  This may be useful to monitor convergence
dets <- apply(res$theta$a.cov, 3, det)
dets.df <- data.frame(x=1:length(dets), dets)

dev.new(height = 5, width = 6)
ggplot(data=dets.df, aes(x=x, y=dets)) +
  geom_point(shape=16, size=4, color="blue", alpha=0.8) +
  labs(x="Iteration", y="Determinant") +
  theme_grey(base_size = 17) +
  scale_x_continuous(breaks=1:length(dets)) +
  scale_y_continuous(breaks=pretty(dets)) 


############## Summary stats from AMIS-INLA

####### Summary stats (mean and variance matrix of IS simulations)
res$theta$a.mu[nrow(res$theta$a.mu), ]
res$theta$a.cov[, , nrow(res$theta$a.mu)]


####### Summary stats using inla.zmarginal

estimations <- lapply(res$margs, inla.zmarginal)
beta_est <- as.vector(estimations[[1]])
precd_est <- as.vector(estimations[[2]])
gamma_est <- as.vector(estimations[[3]])
precu_est <- as.vector(estimations[[4]])


summ <- rbind(beta_est[c(1,2,3,7)],precd_est[c(1,2,3,7)], gamma_est[c(1,2,3,7)], precu_est[c(1,2,3,7)])
summ.is <- apply(unname(summ),c(1,2), FUN=as.numeric)
rownames(summ.is) <- c("beta", "precd", "gamma", "precu")
colnames(summ.is) <- c("Mean", "SD","2.5%","97.5%")

summ.is


############## Summary stats from MCMC

# Load coda samples
#load("coda_sleepstudy_rcoef.Rda")
source("../bugs_models/to_coda_samples.R")
coda_sleepstudy_rcoef <- to_coda_samples(js)

# Show dim. names
dimnames(coda_sleepstudy_rcoef[[1]])

# Index to extract data
idx_vars <- as.integer(sapply(rownames(summ.is), function(X) {
  which(X == dimnames(coda_sleepstudy_rcoef[[1]])[[2]])
}))

# Traceplots and density plots for the estimated parameters

dev.new(height = 6, width = 9)
traplot(coda_sleepstudy_rcoef, parms=c("beta","precd","gamma", "precu"))
dev.new(height = 6, width = 9)
denplot(coda_sleepstudy_rcoef, parms=c("beta","precd","gamma", "precu"))

# Show estimated parameters (Mean, SD and 95% Credible Interval)
summ.mcmc <- cbind(summary(coda_sleepstudy_rcoef)$statistics[idx_vars, c(1,2)],summary(coda_sleepstudy_rcoef)$quantiles[idx_vars, c(1,5)])
summ.mcmc


############## Summary stats from AMIS-INLA and MCMC methods (Mean and SD)
cbind(summ.is[,c(1,2)],summ.mcmc[,c(1,2)])


############## Plot posterior densities

#### Parameter beta

den.mcmc <- density(js$beta[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[1]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p1 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(beta[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 


#### Parameter precd

den.mcmc <- density(js$precd[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[2]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p2 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(tau[beta]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 


#### Parameter gamma
den.mcmc <- density(js$gamma[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[3]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p3 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(gamma), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 

#### Parameter prec_u
den.mcmc <- density(js$precu[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[4]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p4 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(tau[u]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],15)) +
  theme_grey(base_size = 15) 


dev.new(height = 10, width = 10)
ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
