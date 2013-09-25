	
splitEMInit <- function(y, g, conds, lib.size,
	lib.type, alg.type, fixed.lambda, equal.proportions, 
	prev.labels, prev.probaPost, init.runs, init.iter, verbose) {

	## g is the new number of clusters IN ADDITION TO FIXED LAMBDA
	## NB: This function inspiried by init2.k() function from poisson.glm.mix package
	## (written by Panos Papastamoulis)

	n <- dim(y)[1]
	unique.labels <- unique(prev.labels)
	## Check whether any of the clusters has one or zero observations & remove from consideration
	tab <- table(prev.labels)
	if(length(which(tab < 2)) > 0) {
		unique.labels <- names(tab)
	}
	
	K <- g
	if(class(fixed.lambda) == "list") {
		K <- g + length(fixed.lambda);	
	}


	if(init.runs < 1) {
		cat("Exhaustive search is implemented.\n")
		m <- length(unique(prev.labels))
	}
	if(init.runs > 0) m <- init.runs;
#	if(init.runs >= K) m <- K - 1;
	psim <- matrix(NA, nrow = m, ncol = K)
	prev.K <- K - 1
	
	LL.all <- rep(NA, m)
	init.all <- vector("list", m)

	for(iter in 1:m) {
		if(verbose == TRUE) cat("Split small-em run:", iter, "\n");

		## Randomly choose a cluster to split
		if(init.runs < 1) {
			cluster.choose <- unique(prev.labels)[iter]
			index1 <- which(prev.labels == cluster.choose)
		}
		if(init.runs > 1) {
			cluster.choose <- floor(runif(1)*(K-1)+1)
			index1 <- which(prev.labels == cluster.choose)
		}

		## Random selection of observations within splitted component
		u.numbers <- runif(length(index1))
		t <- matrix(0, nrow = n, ncol = K)
		t[,1:prev.K] <- prev.probaPost
		t[index1,g] <- t[index1,cluster.choose] * u.numbers
		t[index1,cluster.choose] <- t[index1,cluster.choose] * (1-u.numbers)

		## Smoothing t values
		epsilon <- 1e-10
    		maxcut <- 1 - epsilon; mincut <- epsilon
   		t <- apply(t, 2, pmax, mincut)
   		t <- apply(t, 2, pmin, maxcut)
		t <- t/rowSums(t)

		## Initialize EM-algorithm with new z values
		initialize <- probaPostInit(y = y, g = g, lib.size = lib.size,
			lib.type = lib.type, conds = conds, alg.type = alg.type,
			fixed.lambda = fixed.lambda,
			equal.proportions = equal.proportions, probaPost.init = t,
			init.iter = init.iter, verbose = verbose)
		LL.all[iter] <- initialize$log.like
		init.all[[iter]] <- initialize
	}
	
	init.select <- which(LL.all == max(LL.all, na.rm = TRUE))[1]
	final.init <- init.all[[init.select]]
	lambda.init <- final.init$lambda
	pi.init <- final.init$pi

	return(list(lambda.init = lambda.init, pi.init = pi.init))
}

