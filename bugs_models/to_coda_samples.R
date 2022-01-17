# jags.samples() -> coda.samples()
# BAsed on rjags::coda.samples()
#
# out: Output from jags.samples()

# function (model, variable.names = NULL, n.iter, thin = 1, na.rm = TRUE, 
to_coda_samples <- function (out) #, variable.names = NULL, n.iter, thin = 1, na.rm = TRUE, 
    #...) 
{
    start <- 1 #model$iter() + thin
    variables.names <- names(out)
    n.iter <- as.integer(dim(js[[1]])[2])
    thin <- 1
    na.rm <- TRUE
    n_chains <- as.integer(dim(js[[1]])[3]) 
    #out <- jags.samples(model, variable.names, n.iter, thin, 
    #    type = "trace", ...)
    ans <- vector("list", n_chains) #vector("list", nchain(model))
    for (ch in 1:n_chains) {
        ans.ch <- vector("list", length(out))
        vnames.ch <- NULL
        for (i in seq(along = out)) {
            varname <- names(out)[[i]]
            d <- dim(out[[i]])
            if (length(d) < 3) {
                stop("Invalid dimensions for sampled output")
            }
            vardim <- d[1:(length(d) - 2)]
            nvar <- prod(vardim)
            niter <- d[length(d) - 1]
            nchain <- d[length(d)]
            values <- as.vector(out[[i]])
            var.i <- matrix(NA, nrow = niter, ncol = nvar)
            for (j in 1:nvar) {
                var.i[, j] <- values[j + (0:(niter - 1)) * nvar + 
                  (ch - 1) * niter * nvar]
            }
            vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
            ans.ch[[i]] <- var.i
        }
        ans.ch <- do.call("cbind", ans.ch)
        colnames(ans.ch) <- vnames.ch
        ans[[ch]] <- mcmc(ans.ch, start = start, thin = thin)
    }
    if (isTRUE(na.rm)) {
        all.missing <- sapply(ans, function(x) {
            apply(is.na(x), 2, any)
        })
        drop.vars <- if (is.matrix(all.missing)) {
            apply(all.missing, 1, any)
        }
        else {
            any(all.missing)
        }
        ans <- lapply(ans, function(x) return(x[, !drop.vars, 
            drop = FALSE]))
    }
    mcmc.list(ans)
}

