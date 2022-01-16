# Simulate data
n <- 500

#Set seed
set.seed(1) 

# Covariates
x_mean1 <- runif(n, 10, 20)
x_prec1 <- runif(n, 0, 20)

# Standardize covariate
x_prec1 <- (x_prec1 - mean(x_prec1))/ sd(x_prec1)


# Compute precisions
# gamma_0 = 0, gamma_1 = 0.5
#size <- exp(0 + 0.5 * x_prec1)
size <- exp(0 + 5 * x_prec1)
summary(size)

# Compute predlin
# beta_0 = 1, beta_1 = 0.25
predlin <- 1 + 0.25 * x_mean1 
summary(predlin)
y <- rnbinom(n, size = size, mu = exp(predlin))

# Summary of simulated data
summary(y)
hist(y)

# Create data.frame with data
d <- data.frame(y, x_mean1, x_prec1,
  idx = 1:n #Index for random effect
)

save(file = "data.RData", list = ls())
