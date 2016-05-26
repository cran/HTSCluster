logLikePoisMix <- 
function (y, mean, pi) 
{
    if (length(mean) != length(pi)) 
        stop(paste(sQuote("mean"), "must be a list of the same length as", 
            sQuote("pi")))
    g <- length(pi)
    n <- dim(y)[1]
    cols <- dim(y)[2]
    nas <- 0
    y <- matrix(y, nrow = n, ncol = cols)
    logLike <- rep(0, n)
    index <- 1:g
    epsilon <- exp(-720)
    thresh <- -745
    nn <- 0
    logpi <- log(pi)
    ef <- matrix(logpi, nrow = n, ncol = g, byrow = T)
    for (k in 1:g) {
        ef[, k] <- ef[, k] + rowSums(dpois(y, mean[[k]], log = T))
    }
    efmax <- apply(ef, 1, max)
    ef <- ef - efmax
    logLike <- efmax + log(rowSums(exp(ef)))
    return(list(ll = sum(logLike, na.rm = TRUE)))
}
