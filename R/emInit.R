
emInit <- function(y, g, conds, norm, alg.type = "EM", 
	init.runs, init.iter, fixed.lambda, equal.proportions, 
	verbose) {
	if(alg.type != "EM" & alg.type != "CEM")
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(length(alg.type) > 1)
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(init.runs < 1 | length(init.runs) > 1 | round(init.runs) != init.runs) 
		stop(paste(sQuote("init.runs"), "must be a positive integer"))
	if(is.logical(verbose) == FALSE)
		stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
	init.type1 <- "kmeans"
	lambda.init.all <- vector("list", init.runs)
	pi.init.all <- vector("list", init.runs)
	criterion.all <- rep(NA, init.runs)
	for(start in 1:init.runs) {
		em.init <- PoisMixClus(y = y, g = g, norm=norm, conds = conds, 
			init.type = init.type1, alg.type = alg.type, iter = init.iter,
			fixed.lambda = fixed.lambda, equal.proportions = equal.proportions,
			wrapper=FALSE)
		lambda.init.all[[start]] <- em.init$lambda
		pi.init.all[[start]] <- em.init$pi
		criterion.all[start] <- em.init$log.like
		if(verbose == TRUE) print(paste("Initialization:", start))
	}
	## If all criterion values are equal to NaN, then arbitrarily choose the first one
	if(sum(is.na(criterion.all)) == length(criterion.all)) {
		final.choice <- 1;
	}
	## If two of the criterion values are exactly the same, pick only the first
	if(sum(is.na(criterion.all)) != length(criterion.all)) {
		final.choice <- which(criterion.all == min(criterion.all, na.rm = TRUE))[1]
	}
	lambda.init <- lambda.init.all[[final.choice]]
	pi.init <- pi.init.all[[final.choice]]
	return(list(pi.init = pi.init, lambda.init = lambda.init))
}
