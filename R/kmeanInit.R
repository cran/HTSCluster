kmeanInit <- function(y, g, conds, norm, fixed.lambda,
	equal.proportions) {

	n <- dim(y)[1];cols <- dim(y)[2];
	y <- as.matrix(y, nrow = n, ncol = cols)
	d <- length(unique(conds))
	r <- as.vector(table(conds))
	w <- rowSums(y)
	
	## Only calculate s values if they are not provided
	if(length(norm) == 1) {
		if(norm == "none") 	s <- rep(1, cols);
		if(norm == "TC") s <- colSums(y) / sum(as.numeric(y));
		if(norm == "UQ") s <- apply(y, 2, quantile, 0.75) / sum(apply(y, 2, quantile, 0.75));
		if(norm == "Med") s <- apply(y, 2, median) / sum(apply(y, 2, median));
		if(norm == "DESeq") {
			## Code from DESeq, v1.8.3
			loggeomeans <- rowMeans(log(y))
			s <- apply(y, 2, function(x) exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
			s <- s / sum(s)
		}
		if(norm == "TMM") {
			f <- calcNormFactors(as.matrix(y), method = "TMM")
			s <- colSums(y)*f / sum(colSums(y)*f)
		} 
	}
	if(length(norm) == length(conds)) {
		s <- norm / sum(norm)
	}
	s.dot <- rep(NA, d) 
	for(j in 1:d) {
		s.dot[j] <- sum(s[which(conds == (unique(conds))[j])])
	}

	## Use K-means to create initial partition
	g.init <- kmeans(y / w, g)
	partition <- g.init$cluster
	partition.mat <- matrix(0, nrow = n, ncol = g)
	for(i in 1:n) partition.mat[i,partition[i]] <- 1;

	## Calculate lambda.init and p.init
	denom <- colSums(partition.mat * w)
	lambda.init <- matrix(NA, nrow = d, ncol = g)
	for(j in 1:d) {
		denom.bis <- denom * s.dot[j]
		num <- colSums(partition.mat *
			rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
		lambda.init[j,] <- num / denom.bis
	}
	if(inherits(fixed.lambda, "list")) {
		index <- lapply(fixed.lambda, function(x) {
			colSums((lambda.init - x)^2)})
		index.order <- lapply(index, order)
		fixed <- c()
		## Pick which components correspond to fixed components
		## Assign fixed lambdas to "closest" component (Euclidean distance),
		## starting with the first fixed lambda
		for(ll in 1:length(fixed.lambda)) {
			while(index.order[[ll]][1] %in% fixed) {
				index.order[[ll]] <- index.order[[ll]][-1]
			}
			lambda.init[,index.order[[ll]][1]] <- fixed.lambda[[ll]]	
			fixed <- c(fixed, index.order[[ll]][1])
		}
		## Rearrange lambda so that fixed values are in the first columns
		lambda.init <- cbind(lambda.init[,fixed], lambda.init[,-fixed])
	}
	if(equal.proportions == FALSE) {
		pi.init <- as.vector(table(g.init$cluster)/n)
		## Rearrange pi so that it corresponds to lambda
		if(is.list(fixed.lambda) == TRUE) {
			pi.init <- c(pi.init[fixed], pi.init[-fixed])
		}
	}
	if(equal.proportions == TRUE) {
		pi.init <- rep(1/g, g)
	}
	return(list(pi.init = pi.init, lambda.init = lambda.init))
}

