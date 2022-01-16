# Gaussian model with hierarchy on the precisions

# Simulate data
# Number of groups
n_groups <- 5
# Number of observations per group
ni <- rep(500, n_groups)
# Number of data
n <- sum(ni)

#Set seed
set.seed(1) 

# Covariates
x_prec1 <- runif(n_groups, -1, 1)

# Compute precisions
# Random effects for precisions
g_reff <- rnorm(n_groups, sd = 1)
# gamma_0 = 0, gamma_1 = 0.5
g_precs <- exp(0 + 5 * x_prec1 +  g_reff)
summary(g_precs)

# Fixed effects versus fixed eff. + error
plot(exp(0 + 5 * x_prec1), g_precs)

# Compute predlin
x_mean1 <- runif(n)
predlin <- 1 + 0.25 * x_mean1
summary(predlin)

# Sample response
y <- rnorm(n, predlin, sd = 1 / sqrt(rep(g_precs, ni)))

# Summary of simulated values
summary(y)
hist(y)

# Create data.frame with data
d <- data.frame(y, x_mean1, g_x_prec1 = rep(x_prec1, ni),
  grp = rep(1:n_groups, ni), # Group index
  idx = 1:n #Individual index
)


# Histogram of within group variances
pdf(file = "hist.pdf")
hist(as.vector(with(d, by(y, grp, var))))
dev.off()

# Boxplot
pdf(file = "boxplot.pdf")
boxplot(y ~ grp, data = d)
dev.off()

save(file = "data.RData", list = ls())
