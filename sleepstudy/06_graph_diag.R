#### Graphical diagnostics

####### Load results
load("AMIS.RData")

####### Re-scaled weights
ww <- res$weight / sum(res$weight)
summary(ww)


# Diagnostic plot
# ww: weights
# smp: Sample (variables by column)

### Part 1
d_plot <- function(ww, smp ) {
  
  # Number of variables
  n.var <- ncol(smp)
  n.samples <- nrow(smp)
  
  for(i in 1:(n.var/2)) {
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

pdf("graph_diag_sleepstudy1.pdf")
par(mfrow = c(3,3))
d_plot(ww, res$eta)
dev.off()



### Part 2
d_plot <- function(ww, smp ) {
  
  # Number of variables
  n.var <- ncol(smp)
  n.samples <- nrow(smp)
  
  for(i in (n.var/2 +1):n.var) {
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


pdf("graph_diag_sleepstudy2.pdf")
par(mfrow = c(3,3))
d_plot(ww, res$eta)
dev.off()
