# Create data objects
#Load libraries
library("spdep")

# Load original data
# This is NOT PROVIDED due to confidentiality issues
load("../Colombia/colombia.Rda")
#head(colombia@data)

# Simulate number of deaths to replace the original ones
set.seed(1234)
sim_data <- rmultinom(1, sum(colombia$Ndeaths), colombia$Ndeaths / sum(colombia$Ndeaths))[, 1]


# Compare to original data
table(sim_data - colombia$Ndeaths)

# Replace original data with simulated data
colombia$Ndeaths <- sim_data

### Order 1 Contiguity spatial weights and spatial lags
nb.list <- poly2nb(colombia, queen=TRUE) # Construct neighbors list
w.list <- nb2listw(nb.list, style="W", zero.policy=TRUE)
colombia$rates <- colombia$Ndeaths/colombia$N0years*1000
colombia$lag.rates <- lag.listw(w.list,colombia$rates,zero.policy=TRUE)


n <- nrow(colombia)
y <- colombia$Ndeaths
offset1 <- colombia$N0years
x_mean1 <- colombia$lag.rates
x_prec1 <- colombia$NBI
d <- data.frame(y, offset1, x_mean1, x_prec1, idx = 1:n) 



save(file = "data.RData", list = ls())
