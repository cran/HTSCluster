logLikePoisMix <-
function(y, mean, pi) {

	if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
		stop(paste(sQuote("y"), "must be a matrix"))
	if(min(y) < 0 | sum(round(y)) != sum(y)) 
		stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
	if(min(rowSums(y)) == 0)
		stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
	if(length(mean) != length(pi)) 
		stop(paste(sQuote("mean"), "must be a list of the same length as", sQuote("pi")))

	g <- length(pi)
	n <- dim(y)[1]; cols <- dim(y)[2];
	y <- matrix(y, nrow = n, ncol = cols)

	## Using calculation trick from Panos
	h <- matrix(NA, nrow = n, ncol = g)
	for(k in 1:g) {
		h[,k] <- log(pi[k]) + rowSums(log((dpois(y, mean[[k]]))))
	}
	h.star <- apply(h, 1, max)
	logLike <- sum(h.star + log(rowSums(exp(h - h.star))), na.rm = TRUE)

#	logLike <- rep(0, n)
#	for(i in 1:n) {
#		tmp <- 0
#		for(k in 1:g) {
#			dpois.tmp <- dpois(y[i,], mean[[k]][i,])
#			if(smoothing == TRUE) {
#				## Smoothing (0's set to 1e-10, 1's set to 1-1e-10)
#				epsilon <- 1e-10
#				maxcut <- 1-epsilon; mincut <- epsilon;
#				dpois.tmp <- pmax(dpois.tmp, mincut)
#				dpois.tmp <- pmin(dpois.tmp, maxcut);
#			}
#			dpois.tmp <- prod(dpois.tmp)
#			tmp <- tmp + pi[k] * prod(dpois.tmp)
#		}
#		tmp <- ifelse(tmp == 0, NA, tmp)
#		logLike[i] <- log(as.numeric(tmp))
#	}
#	return(list(ll = sum(logLike, na.rm = TRUE)))

	return(list(ll = logLike))
}

