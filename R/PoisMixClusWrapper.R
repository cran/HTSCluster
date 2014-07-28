PoisMixClusWrapper <- function(y, gmin = 1, gmax, conds, lib.size = TRUE, lib.type = "TMM", 
	gmin.init.type = "small-em", init.runs = 1, init.iter = 10, split.init = TRUE, 
	alg.type = "EM", cutoff = 10e-6, iter = 1000, 
	fixed.lambda = NA, equal.proportions = FALSE, verbose = FALSE,
	interpretation = "sum", EM.verbose = FALSE) 
{
	all.results <- vector("list", length = gmax - gmin + 1)
	names(all.results) <- paste("g=", seq(gmin,gmax, 1), sep = "")

	## For gmin, run PoisMixClus with regular small-EM initialization
	cat("Running g =", gmin, "...\n")
	all.results[[1]] <- PoisMixClus(y = y, g = gmin, lib.size = lib.size, 
		lib.type = lib.type, conds = conds, 
		init.type = gmin.init.type, init.runs = init.runs, init.iter = init.iter,
		alg.type = alg.type, cutoff = cutoff, iter = iter, 
		fixed.lambda = fixed.lambda, equal.proportions = equal.proportions, 
		prev.labels = NA, prev.probaPost = NA, verbose = verbose,
 interpretation = interpretation,
		EM.verbose = EM.verbose)

	## For g > gmin, run PoisMixClus with Panos-like init using previous results
	index <- 2
	if(gmax > gmin) {
		if(split.init == TRUE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				prev.labels <- all.results[[K-1]]$labels
				prev.probaPost <- all.results[[K-1]]$probaPost
				all.results[[index]] <- PoisMixClus(y = y, g = K, lib.size = lib.size, 
					lib.type = lib.type, conds = conds, 
					init.type = "split.small-em", 
					alg.type = alg.type, cutoff = cutoff, iter = iter, 
					fixed.lambda = fixed.lambda, 
					equal.proportions = equal.proportions, 
					prev.labels = prev.labels, prev.probaPost = prev.probaPost,
					init.runs = init.runs, init.iter = init.iter, verbose = verbose,
										interpretation = interpretation, EM.verbose = EM.verbose)
				index <- index + 1
			}
		}
		if(split.init == FALSE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				all.results[[index]] <- PoisMixClus(y = y, g = K, lib.size = lib.size, 
					lib.type = lib.type, conds = conds, 
					init.type = gmin.init.type, 
					alg.type = alg.type, cutoff = cutoff, iter = iter, 
					fixed.lambda = fixed.lambda, 
					equal.proportions = equal.proportions, 
					prev.labels = NA, prev.probaPost = NA, init.runs = init.runs,
					init.iter = init.iter, verbose = verbose, interpretation = interpretation,
					EM.verbose = EM.verbose)
				index <- index + 1
			}
		}
	}

	logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
	ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
	ICL.choose <- which(ICL.all == max(ICL.all, na.rm = TRUE))
	select.results <- all.results[[ICL.choose]]
	RESULTS <- list(logLike.all = logLike.all, ICL.all = ICL.all,
		select.results = select.results, all.results = all.results)
	class(RESULTS) <- "HTSClusterWrapper"
	return(RESULTS)
}
