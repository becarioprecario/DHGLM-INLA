############################### Summary of results 
############################### Poisson model (Colombia data)

library(INLA)
library(mcmcplots)
library(ggplot2)
library(coda)
library(rjags)
library(ggpubr) 


####### Load results
load("AMIS_PO.RData")
load("JAGS_PO.RData")


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
  scale_x_continuous(breaks=1:length(dets))


############## Summary stats from AMIS-INLA

####### Summary stats (mean and variance matrix of IS simulations)
res$theta$a.mu[nrow(res$theta$a.mu), ]
res$theta$a.cov[, , nrow(res$theta$a.mu)]


####### Summary stats using inla.zmarginal

# beta
beta_est <- lapply(res$margs, inla.zmarginal)
beta1_est <- as.vector(beta_est[[1]])
rho_est <- as.vector(beta_est[[2]])

# gamma
gamma_est <- lapply(res$eta_kern, inla.zmarginal)
gamma1_est <- as.vector(gamma_est[[1]])
gamma2_est <- as.vector(gamma_est[[2]])

summ <- rbind(beta1_est[c(1,2,3,7)],rho_est[c(1,2,3,7)], gamma1_est[c(1,2,3,7)], gamma2_est[c(1,2,3,7)])
summ.is <- apply(unname(summ),c(1,2), FUN=as.numeric)
rownames(summ.is) <- c("beta", "rho", "gamma0", "gamma1")
colnames(summ.is) <- c("Mean", "SD","2.5%","97.5%")

summ.is



############## Summary stats from MCMC

# Load coda samples
#load("coda_colpoi.Rda")
source("../bugs_models/to_coda_samples.R")
coda_colpoi <- to_coda_samples(js) 

# Traceplots and density plots for the estimated parameters
dev.new(height = 6, width = 9)
traplot(coda_colpoi, parms=c("beta","gamma"))
dev.new(height = 6, width = 9)
denplot(coda_colpoi, parms=c("beta","gamma"))

# Show estimated parameters (Mean, SD and 95% Credible Interval)
summ.mcmc <- cbind(summary(coda_colpoi)$statistics[c(1:4),c(1,2)],summary(coda_colpoi)$quantiles[c(1:4),c(1,5)])
rownames(summ.mcmc) <- c("beta", "rho", "gamma0", "gamma1")
summ.mcmc


############## Summary stats from AMIS-INLA and MCMC methods (Mean and SD)
cbind(summ.is[,c(1,2)],summ.mcmc[,c(1,2)])


############## Plot posterior densities

#### Parameter beta_0 
den.mcmc <- density(js$beta[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[1]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p1 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(beta), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 


#### Parameter rho

den.mcmc <- density(js$beta[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$margs[[2]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p2 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(rho), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 

#### Parameter gamma_0  
den.mcmc <- density(js$gamma[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$eta_kern[[1]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p3 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(gamma[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 


#### Parameter gamma_1
den.mcmc <- density(js$gamma[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")
den.is <- res$eta_kern[[2]]

join_den <- rbind(den.mc,den.is)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA", nrow(den.is)))

p4 <- ggplot(data=join_den, aes(x=x, y=y, fill=dataset)) + 
  geom_area(alpha=0.5) +
  scale_fill_manual(name="Algorithm", values=c("#DF1312", "#131597")) + 
  labs(x=expression(gamma[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  theme_grey(base_size = 15) 

dev.new(height = 10, width = 10)
ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
