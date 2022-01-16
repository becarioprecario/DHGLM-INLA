# Simulate data
n <- 1000

#Set seed
set.seed(1) 

# Covariates
x_mean1 <- runif(n)
x_prec1 <- rnorm(n)

# Compute precisions
# gamma_0 = 0, gamma_1 = 0.5
precs <- exp(0 + 0.5 * x_prec1)
summary(precs)

# Compute predlin
reff <- rnorm(n, mean = 0, sd = sqrt(1 / precs))
# beta_0 = 1, beta_1 = 0.25
predlin <- 1 + 0.25 * x_mean1 + reff
summary(predlin)
y <- rpois(n, exp(predlin))

# Summary of simulated data
summary(y)
hist(y)

# Create data.frame with data
d <- data.frame(y, x_mean1, x_prec1,
  idx = 1:n #Index for random effect
)


save(file = "data.RData", list = ls())

