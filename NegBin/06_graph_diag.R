#### Graphical diagnostics

####### Load results
load("AMIS.RData")

####### Re-scaled weights
ww <- res$weight / sum(res$weight)
summary(ww)


# Diagnostic plot
# ww: weights
# smp: Sample (variables by column)
d_plot <- function(ww, smp ) {
  
  # Number of variables
  n.var <- ncol(smp)
  n.samples <- nrow(smp)
  
  for(i in 1:n.var) {
    # Reorder weights
    idx <- order(smp[, i])
    plot(1:n.samples / n.samples, cumsum(ww[idx]),
         main = bquote(gamma[~ .(i-1)]), type = "l",
         xlab = "Theoretical cumulative distribution",
         ylab = "Empirical cumulative distribution")
    abline(0, 1, lty = 2)
  }
}


pdf("graph_diag_simnegbin.pdf", width=7, height=4)
par(mfrow = c(1, 2))
d_plot(ww, res$eta)
dev.off()
