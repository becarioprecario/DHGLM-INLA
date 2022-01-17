############################### Summary of results 
############################### Gaussian (simulated data)
# True values: beta_0 = 1, beta_1 = 0.25, gamma_0 = 0, gamma_1 = 5, precu = 1


library(INLA)
library(mcmcplots)
library(ggplot2)
library(coda)
library(rjags)
library(ggpubr) 


####### Load results
load("JAGS.RData")

amis <- lapply(paste("AMIS", c("", "B", "C", "D", "B10","B55K"), ".RData", sep = ""),
               function(X) {mget(load(X))})
# Transform marginals 'prec internal scale' -> 'prec.'
for(i in 1:length(amis)) {
   amis[[i]]$res$margs[[5]]<- data.frame(inla.tmarginal(exp, amis[[i]]$res$margs[[6]]))
}

####### Computing times
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



########### Determinant of variance matrix of importance distribution
#  This may be useful to monitor convergence
# There are actual zeros for cases A and C
dev.new(height = 5, width = 6)
plot(log(apply(amis[[1]]$res$theta$a.cov, 3, det)), type = "n")
lapply(1:6, function(X) {
  lines(log(apply(amis[[X]]$res$theta$a.cov, 3, det)), lty = X)
})


############## Summary stats from AMIS-INLA

####### Summary stats (mean and variance matrix of IS simulations)
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

####### Summary stats using inla.zmarginal

estimations <- lapply(amis, function(X) {
  with(X, lapply(res$margs, inla.zmarginal))
})

mult_amis <- function(X){
  ## amis1
  estimations <- X
  
  # beta
  beta1_est <- unlist(estimations[[1]])
  beta2_est <- unlist(estimations[[2]])
  
  # gamma
  gamma1_est <- unlist(estimations[[3]])
  gamma2_est <- unlist(estimations[[4]])
  
  #precu
  precu_est <- unlist(estimations[[5]])
  
  summ <- rbind(beta1_est[c(1,2,3,7)],beta2_est[c(1,2,3,7)], gamma1_est[c(1,2,3,7)], gamma2_est[c(1,2,3,7)],precu_est[c(1,2,3,7)])
  summ.is <- apply(unname(summ),c(1,2), FUN=as.numeric)
  rownames(summ.is) <- c("beta0", "beta1", "gamma0", "gamma1","precu")
  colnames(summ.is) <- c("Mean", "SD", "95 CI","0")
  
  return(summ.is)
}


summ.is1 <- mult_amis(estimations[[1]])
summ.is1

summ.is2 <- mult_amis(estimations[[2]])
summ.is2

summ.is3 <- mult_amis(estimations[[3]])
summ.is3

summ.is4 <- mult_amis(estimations[[4]])
summ.is4

summ.is5 <- mult_amis(estimations[[5]])
summ.is5

summ.is6 <- mult_amis(estimations[[6]])
summ.is6


############## Summary stats from MCMC

# Load coda samples
#load("coda_simgauss5.Rda")
source("../bugs_models/to_coda_samples.R")
coda_simgauss5 <- to_coda_samples(js) #[c("beta", "gamma", "precu")])


# Traceplots and density plots for the estimated parameters
dev.new(height = 6, width = 9)
traplot(coda_simgauss5)
dev.new(height = 6, width = 9)
denplot(coda_simgauss5)

# Show estimated parameters (Mean, SD and 95% Credible Interval)
summ.mcmc <- cbind(summary(coda_simgauss5)$statistics[c(1:4,10),c(1,2)],summary(coda_simgauss5)$quantiles[c(1:4,10),c(1,5)])
summ.mcmc


############## Plot posterior densities

##### All scenarios

#### Parameter beta_0 

den.mcmc <- density(js$beta[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[1]]
den.is2 <- amis[[2]]$res$margs[[1]]
den.is3 <- amis[[3]]$res$margs[[1]]
den.is4 <- amis[[4]]$res$margs[[1]]
den.is5 <- amis[[5]]$res$margs[[1]]
den.is6 <- amis[[6]]$res$margs[[1]]


join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4,den.is5,den.is6)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)),
                      rep("AMIS-INLA5", nrow(den.is5)),rep("AMIS-INLA6", nrow(den.is6)))

p1 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c("#32A251","#835B82","#FC7715","#B71E42","#FFC000", "#686767","#131597")) + 
  labs(x=expression(beta[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 

#### Parameter beta_1

den.mcmc <- density(js$beta[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[2]]
den.is2 <- amis[[2]]$res$margs[[2]]
den.is3 <- amis[[3]]$res$margs[[2]]
den.is4 <- amis[[4]]$res$margs[[2]]
den.is5 <- amis[[5]]$res$margs[[2]]
den.is6 <- amis[[6]]$res$margs[[2]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4,den.is5,den.is6)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)),
                      rep("AMIS-INLA5", nrow(den.is5)),rep("AMIS-INLA6", nrow(den.is6)))

#dev.new(height = 5, width = 7)
p2 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c("#32A251","#835B82","#FC7715","#B71E42","#FFC000", "#686767","#131597")) + 
  labs(x=expression(beta[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0.25, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_0  
den.mcmc <- density(js$gamma[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[3]]
den.is2 <- amis[[2]]$res$margs[[3]]
den.is3 <- amis[[3]]$res$margs[[3]]
den.is4 <- amis[[4]]$res$margs[[3]]
den.is5 <- amis[[5]]$res$margs[[3]]
den.is6 <- amis[[6]]$res$margs[[3]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4,den.is5,den.is6)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)),
                      rep("AMIS-INLA5", nrow(den.is5)),rep("AMIS-INLA6", nrow(den.is6)))

#dev.new(height = 5, width = 7)
p3 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c("#32A251","#835B82","#FC7715","#B71E42","#FFC000", "#686767","#131597")) + 
  labs(x=expression(gamma[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_1
den.mcmc <- density(js$gamma[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[4]]
den.is2 <- amis[[2]]$res$margs[[4]]
den.is3 <- amis[[3]]$res$margs[[4]]
den.is4 <- amis[[4]]$res$margs[[4]]
den.is5 <- amis[[5]]$res$margs[[4]]
den.is6 <- amis[[6]]$res$margs[[4]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4,den.is5,den.is6)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)),
                      rep("AMIS-INLA5", nrow(den.is5)),rep("AMIS-INLA6", nrow(den.is6)))

#dev.new(height = 5, width = 7)
p4 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c("#32A251","#835B82","#FC7715","#B71E42","#FFC000", "#686767","#131597")) + 
  labs(x=expression(gamma[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=5, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter prec_u
den.mcmc <- density(js$precu[, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[5]]
den.is2 <- amis[[2]]$res$margs[[5]]
den.is3 <- amis[[3]]$res$margs[[5]]
den.is4 <- amis[[4]]$res$margs[[5]]
den.is5 <- amis[[5]]$res$margs[[5]]
den.is6 <- amis[[6]]$res$margs[[5]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4,den.is5,den.is6)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)),
                      rep("AMIS-INLA5", nrow(den.is5)),rep("AMIS-INLA6", nrow(den.is6)))

p5 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c("#32A251","#835B82","#FC7715","#B71E42","#FFC000", "#686767","#131597")) + 
  labs(x=expression(tau[u]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1]-5,range(den.mc$x)[2])+5) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 


dev.new(height = 13, width = 11)
ggarrange(p1,p2,p3,p4,p5, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")


##### Scenarios 1 to 4

#### Parameter beta_0 

den.mcmc <- density(js$beta[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[1]]
den.is2 <- amis[[2]]$res$margs[[1]]
den.is3 <- amis[[3]]$res$margs[[1]]
den.is4 <- amis[[4]]$res$margs[[1]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)))

p1 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",4), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:4,1)) +
  labs(x=expression(beta[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter beta_1

den.mcmc <- density(js$beta[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[2]]
den.is2 <- amis[[2]]$res$margs[[2]]
den.is3 <- amis[[3]]$res$margs[[2]]
den.is4 <- amis[[4]]$res$margs[[2]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)))

#windows(height = 5, width = 7)
p2 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",4), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:4,1)) +
  labs(x=expression(beta[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0.25, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_0  
den.mcmc <- density(js$gamma[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[3]]
den.is2 <- amis[[2]]$res$margs[[3]]
den.is3 <- amis[[3]]$res$margs[[3]]
den.is4 <- amis[[4]]$res$margs[[3]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)))

#windows(height = 5, width = 7)
p3 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",4), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:4,1)) +
  labs(x=expression(gamma[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_1
den.mcmc <- density(js$gamma[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[4]]
den.is2 <- amis[[2]]$res$margs[[4]]
den.is3 <- amis[[3]]$res$margs[[4]]
den.is4 <- amis[[4]]$res$margs[[4]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)))

#windows(height = 5, width = 7)
p4 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",4), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:4,1)) +
  labs(x=expression(gamma[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=5, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter prec_u
den.mcmc <- density(js$precu[, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[1]]$res$margs[[5]]
den.is2 <- amis[[2]]$res$margs[[5]]
den.is3 <- amis[[3]]$res$margs[[5]]
den.is4 <- amis[[4]]$res$margs[[5]]

join_den <- rbind(den.mc,den.is1,den.is2,den.is3,den.is4)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA1", nrow(den.is1)),rep("AMIS-INLA2", nrow(den.is2)),
                      rep("AMIS-INLA3", nrow(den.is3)),rep("AMIS-INLA4", nrow(den.is4)))

p5 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",4), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:4,1)) +
  labs(x=expression(tau[u]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1]-5,range(den.mc$x)[2])+5) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 


dev.new(height = 13, width = 11)
ggarrange(p1,p2,p3,p4,p5, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")



##### Scenarios 5 and 6

#### Parameter beta_0 

den.mcmc <- density(js$beta[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[5]]$res$margs[[1]]
den.is2 <- amis[[6]]$res$margs[[1]]


join_den <- rbind(den.mc,den.is1,den.is2)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA5", nrow(den.is1)),rep("AMIS-INLA6", nrow(den.is2)))

#dev.new(height = 7, width = 9)
p1 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",2), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:2,1)) +
  labs(x=expression(beta[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter beta_1

den.mcmc <- density(js$beta[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[5]]$res$margs[[2]]
den.is2 <- amis[[6]]$res$margs[[2]]


join_den <- rbind(den.mc,den.is1,den.is2)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA5", nrow(den.is1)),rep("AMIS-INLA6", nrow(den.is2)))

#dev.new(height = 5, width = 7)
p2 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",2), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:2,1)) +
  labs(x=expression(beta[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0.25, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_0  
den.mcmc <- density(js$gamma[1, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[5]]$res$margs[[3]]
den.is2 <- amis[[6]]$res$margs[[3]]


join_den <- rbind(den.mc,den.is1,den.is2)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA5", nrow(den.is1)),rep("AMIS-INLA6", nrow(den.is2)))

#dev.new(height = 5, width = 7)
p3 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",2), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:2,1)) +
  labs(x=expression(gamma[0]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=0, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter gamma_1
den.mcmc <- density(js$gamma[2, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[5]]$res$margs[[4]]
den.is2 <- amis[[6]]$res$margs[[4]]


join_den <- rbind(den.mc,den.is1,den.is2)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA5", nrow(den.is1)),rep("AMIS-INLA6", nrow(den.is2)))

#dev.new(height = 5, width = 7)
p4 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",2), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:2,1)) +
  labs(x=expression(gamma[1]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1],range(den.mc$x)[2])) +
  geom_vline(xintercept=5, lty=2, size=1) +
  theme_grey(base_size = 15) 


#### Parameter prec_u
den.mcmc <- density(js$precu[, , 1])
den.mc <- as.data.frame(cbind(den.mcmc$x, den.mcmc$y))
names(den.mc) <- c("x","y")

den.is1 <- amis[[5]]$res$margs[[5]]
den.is2 <- amis[[6]]$res$margs[[5]]


join_den <- rbind(den.mc,den.is1,den.is2)
join_den$dataset <- c(rep("MCMC", nrow(den.mc)), rep("AMIS-INLA5", nrow(den.is1)),rep("AMIS-INLA6", nrow(den.is2)))

#dev.new(height = 5, width = 7)
p5 <- ggplot(data=join_den, aes(x=x, y=y, col=dataset, linetype=dataset)) +
  geom_line( size=1.3) +
  scale_color_manual(name="Algorithm", values=c(rep("#DF1312",2), "#131597")) + 
  scale_linetype_manual(name = "Algorithm", values=c(1:2,1)) +
  labs(x=expression(tau[u]), y="") +
  scale_x_continuous(limits = c(range(den.mc$x)[1]-5,range(den.mc$x)[2])+5) +
  geom_vline(xintercept=1, lty=2, size=1) +
  theme_grey(base_size = 15) 


dev.new(height = 13, width = 11)
ggarrange(p1,p2,p3,p4,p5, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
