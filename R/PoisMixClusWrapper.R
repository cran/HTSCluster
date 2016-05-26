PoisMixClusWrapper <- function(y, gmin=1, gmax, conds, norm="TMM", 
	gmin.init.type = "small-em", init.runs=1, init.iter=10, split.init=TRUE, 
	alg.type="EM", cutoff = 10e-6, iter = 1000, 
	fixed.lambda=NA, equal.proportions = FALSE, verbose = FALSE,
	interpretation = "sum", EM.verbose=FALSE, subset.index = NA) 
{
	if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
		stop(paste(sQuote("y"), "must be a matrix"))
	if(min(y) < 0 | sum(as.numeric(round(y))) != sum(as.numeric(y))) 
		stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
	if(min(rowSums(y)) == 0)
		stop(paste("at least one observation in", sQuote("y"), 
		"contains all 0's and must be removed from the data"))
	if(is.vector(conds) == FALSE | length(conds) != ncol(y))
		stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote("y")))
	if(gmin.init.type != "small-em" & gmin.init.type != "kmeans" & gmin.init.type != "split.small-em") 
		stop(paste(sQuote("gmin.init.type"), "must be one of", dQuote("small-em"), "or", 
			dQuote("kmeans"), "or", dQuote("split.small-em")))
	if(alg.type != "EM" & alg.type != "CEM")
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(length(alg.type) > 1)
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(is.logical(verbose) == FALSE)
		stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
	if(class(fixed.lambda) != "list" & is.na(fixed.lambda[1]) == FALSE)
		stop(paste(sQuote("fixed.lambda"), "must be", dQuote("NA") , "or a list."))
		
	## Grouping columns of y in order of condition (all replicates put together)
	o.ycols <- order(conds)
	y <- y[,o.ycols]
	conds <- conds[o.ycols]
	conds.names <- unique(conds)
	d <- length(unique(conds))
	r <- as.vector(table(conds))
	if(length(rownames(y)) == 0) rn <- 1:nrow(y);
	if(length(rownames(y)) > 0) rn <- rownames(y);
	y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
	rownames(y) <- rn;

	if(length(norm) != 1 & length(norm) != length(conds)) stop(paste(sQuote("norm"), "must be one of
		the following: none, TC, UQ, Med, DESeq, TMM, or a vector of the same length as", sQuote("conds")))
	## If estimated from data, all genes should be used
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
	
	## In case only a subset of data are to be used for analysis
	if(is.na(subset.index)[1] == FALSE) {
		y <- y[subset.index,]
		n <- dim(y)[1];cols <- dim(y)[2]
		w <- rowSums(y)
	}
	if(is.na(subset.index)[1] == TRUE) {
		n <- dim(y)[1];cols <- dim(y)[2]
		w <- rowSums(y)
	}

	all.results <- vector("list", length = gmax - gmin + 1)
	names(all.results) <- paste("g=", seq(gmin,gmax, 1), sep = "")

	## For gmin, run PoisMixClus with regular small-EM initialization
	cat("Running g =", gmin, "...\n")
	all.results[[1]] <- PoisMixClus(y = y, g = gmin, norm = s, conds = conds, 
		init.type = gmin.init.type, init.runs = init.runs, init.iter = init.iter,
		alg.type = alg.type, cutoff = cutoff, iter = iter, 
		fixed.lambda = fixed.lambda, equal.proportions = equal.proportions, 
		prev.labels = NA, prev.probaPost = NA, verbose = verbose,
		interpretation = interpretation,
		EM.verbose = EM.verbose, wrapper=TRUE, subset.index=subset.index)

	## For g > gmin, run PoisMixClus with Panos-like init using previous results
	index <- 2
	if(gmax > gmin) {
		if(split.init == TRUE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				prev.labels <- all.results[[index-1]]$labels
				prev.probaPost <- all.results[[index-1]]$probaPost
				all.results[[index]] <- PoisMixClus(y = y, g = K, norm=s,
					conds = conds, init.type = "split.small-em", 
					alg.type = alg.type, cutoff = cutoff, iter = iter, 
					fixed.lambda = fixed.lambda, 
					equal.proportions = equal.proportions, 
					prev.labels = prev.labels, prev.probaPost = prev.probaPost,
					init.runs = init.runs, init.iter = init.iter, verbose = verbose,
					interpretation = interpretation, EM.verbose = EM.verbose,
					wrapper=TRUE, subset.index=subset.index)
				index <- index + 1
			}
		}
		if(split.init == FALSE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				all.results[[index]] <- PoisMixClus(y = y, g = K,
					norm = s, conds = conds, 
					init.type = gmin.init.type, 
					alg.type = alg.type, cutoff = cutoff, iter = iter, 
					fixed.lambda = fixed.lambda, 
					equal.proportions = equal.proportions, 
					prev.labels = NA, prev.probaPost = NA, init.runs = init.runs,
					init.iter = init.iter, verbose = verbose, interpretation = interpretation,
					EM.verbose = EM.verbose, wrapper=TRUE, subset.index=subset.index)
				index <- index + 1
			}
		}
	}

	logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
	ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
	ICL.choose <- which(ICL.all == max(ICL.all, na.rm = TRUE))
	select.results <- all.results[[ICL.choose]]
	select.results$model.selection <- "ICL"
	
	BIC.all <- unlist(lapply(all.results, function(x) x$BIC))
	BIC.choose <- which(BIC.all == max(BIC.all, na.rm = TRUE))
	select.results2 <- all.results[[BIC.choose]]
	select.results2$model.selection <- "BIC"
	
	# Apply capushe: only if at least 10 models are considered
	if(gmax-gmin+1 <= 10) {
		message("Note: the slope heuristics approach for model selection may only be applied if more than 10 models are fit.")
		DDSE.results <- NA
		Djump.results <- NA
		capushe <- NA
		ResCapushe <- NA
	}
	if(c(gmax-gmin+1) > 10) {
		message("Note: diagnostic plots for results corresponding to model selection via slope heuristics (Djump and DDSE) should be examined to ensure that sufficiently complex models have been considered.")
		Kchoice <- gmin:gmax
		np <- (Kchoice-1) + (length(unique(conds))-1)*(Kchoice)
		mat <- cbind(Kchoice, np/n, np/n, -logLike.all/n)
		ResCapushe <- capushe(mat, n)
		DDSE <- ResCapushe@DDSE@model
		Djump <- ResCapushe@Djump@model
		DDSE.results <- all.results[[paste("g=", DDSE, sep="")]]
		Djump.results <- all.results[[paste("g=", Djump, sep="")]]
		DDSE.results$model.selection <- "DDSE"
		Djump.results$model.selection <- "Djump"
	}
		
	RESULTS <- list(logLike.all = logLike.all, ICL.all = ICL.all,
		capushe = ResCapushe, 
		all.results = all.results,
		DDSE.results = DDSE.results,
		Djump.results = Djump.results,
		ICL.results = select.results,
		BIC.results = select.results2)
	class(RESULTS) <- "HTSClusterWrapper"
	return(RESULTS)
}
