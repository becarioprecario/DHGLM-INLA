# Load and prepare data

library("lme4")
data(sleepstudy)

# Divie times by 1000
sleepstudy$Reaction <- sleepstudy$Reaction / 1000 

# Number of groups
n_groups <- nlevels(sleepstudy$Subject)
# Number of observations per group
ni <- as.vector(table(sleepstudy$Subject))
# Number of data
n <- sum(ni)

# Covariate
x_mean1 <- sleepstudy$Days

# Response
y <- sleepstudy$Reaction
#grp = rep(1:n_groups, ni)
grp <- sleepstudy$Subject

# Create data.frame with data
d <- data.frame(y, x_mean1, 
  grp, # Group index
  idx = 1:n #Individual index
)

save(file = "data.RData", list = ls())
