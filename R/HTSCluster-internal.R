.myloopfxn <- function(k, lambda, w.mat, s.mat, r, n, cols) {
	lambda.mat <- matrix(rep(rep(lambda[,k], times = r), each = n),
		nrow = n, ncol = cols)
	return(w.mat * s.mat * lambda.mat)
}
.myfxn <- function(var1, var2) {
	tmp <- var1 * log(var1/var2) + var2 - var1
	index <- which(var1 == 0)
	tmp[index] <- var2[index]
	rowSums(tmp)
}
.myprobafxn <- function(k, y, pi, mean) {
	pi[k] * exp(rowSums(dpois(y, mean[[k]], log=T)))
}

