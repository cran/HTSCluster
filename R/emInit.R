emInit <-
function(y, g, conds, lib.size, lib.type = "TC", alg.type = "EM", starts = 5, verbose = FALSE) {
if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0)
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(g) != 1)
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(g <= 0 | round(g) != g) 
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(is.vector(conds) == FALSE | length(conds) != ncol(y))
	stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
if(is.logical(lib.size) == FALSE)
	stop(paste(sQuote("libsize"), "must be", dQuote("TRUE"), "(PMM-II) or", 
		dQuote("FALSE"), "(PMM-I)"))
if(lib.type != "TC" & lib.type != "Q" & lib.type != "MedRatio")
	stop(paste(sQuote("lib.type"), "must be one of", dQuote("TC"), "(Total Count),", 
		dQuote("Q"), "(Quantile), or", dQuote("MedRatio"), "(Median Ratio)"))
if(length(lib.type) > 1)
	stop(paste(sQuote("lib.type"), "must be one of", dQuote("TC"), "(Total Count),", 
		dQuote("Q"), "(Quantile), or", dQuote("MedRatio"), "(Median Ratio)"))
if(alg.type != "EM" & alg.type != "CEM")
	stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
if(length(alg.type) > 1)
	stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
if(starts < 1 | length(starts) > 1 | round(starts) != starts) 
	stop(paste(sQuote("starts"), "must be a positive integer"))
if(is.logical(verbose) == FALSE)
	stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
init.type1 <- "kmeans"
lambda.init.all <- vector("list", starts)
pi.init.all <- vector("list", starts)
criterion.all <- rep(NA, starts)
for(start in 1:starts) {
em.init <- PoisMixClus(y, gmin = g, gmax = g, lib.size, lib.type, conds, 
init.type1, alg.type = alg.type, iter = 5)
lambda.init.all[[start]] <- em.init$lambda[[1]]
pi.init.all[[start]] <- em.init$pi[[1]]
criterion.all[start] <- ifelse(alg.type == "EM", em.init$BIC, em.init$ICL)
if(verbose == TRUE) print(paste("Initialization:", start))
}
## If two of the criterion values are exactly the same, pick only the first
final.choice <- which(criterion.all == min(criterion.all))[1]
lambda.init <- lambda.init.all[[final.choice]]
pi.init <- pi.init.all[[final.choice]]
return(list(pi.init = pi.init, lambda.init = lambda.init))
}

