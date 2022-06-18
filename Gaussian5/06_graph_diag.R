#### Graphical diagnostics

library(INLA)

####### Load results
amis <- lapply(paste("AMIS", c("", "B", "C", "D", "B10","B55K"), ".RData", sep = ""),
               function(X) {mget(load(X))})
# Transform marginals 'prec internal scale' -> 'prec.'
for(i in 1:length(amis)) {
  amis[[i]]$res$margs[[5]]<- data.frame(inla.tmarginal(exp, amis[[i]]$res$margs[[6]]))
}

#Re-scaled weights
ww_list <- lapply(amis, function(X) {
  with(X, res$weight / sum(res$weight))
})

# Diagnostic plot
# ww: weights
# smp: Sample (variables by column)

### Part 1
d_plot <- function(ww, smp ) {
  
  # Number of variables
  n.var <- ncol(smp)
  n.samples <- nrow(smp)
  
  for(i in 1:(n.var)) {
    # Reorder weights
    idx <- order(smp[, i])
    plot(1:n.samples / n.samples, cumsum(ww[idx]),
         main = bquote(tau[~ .(i)]), type = "l", 
         cex.lab=1.1, cex.axis=1.1, cex.sub=1.1, cex.main=1.8,
         xlab = "Theoretical cumulative distribution",
         ylab = "Empirical cumulative distribution")
    abline(0, 1, lty = 2)
  }
}

windows(height=9,width=7)
par(mfrow = c(3,2))
d_plot(unlist(ww_list[3]), amis[[3]]$res$eta)

for(i in 1:6){
  path = paste0("graph_diag_gaussian",i,".pdf")
  pdf(path, height=8, width=5)
  par(mfrow = c(3,2))
  d_plot(unlist(ww_list[i]), amis[[i]]$res$eta)
  dev.off()
}


