
`summary.HTSCluster` <-
function (object, cluster.choice = "ICL", ...) 
{
	x <- object
    	if (class(x) != "HTSCluster") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSCluster"), sep = ""), sep = "")
    	}

	if(cluster.choice != "ICL" & cluster.choice != "BIC" & length(cluster.choice) > 1) {
		stop(paste(sQuote("cluster.choice"), sep = ""), " must be one of ",
			paste(dQuote("ICL"), sep = ""), " or ", paste(dQuote("BIC"), sep = ""), 
			" or \n one of the cluster sizes run for object ", paste(sQuote("x"), sep = ""))
	}

	if(cluster.choice == "ICL") {
		probaPost <- x$probaPost.ICL
		g <- x$g.ICL
		labels <- x$labels.ICL
		lambda <- x$lambda.ICL
		pi <- x$pi.ICL
	}
	if(cluster.choice == "BIC") {
		probaPost <- x$probaPost.BIC
		g <- x$g.BIC
		labels <- x$labels.BIC
		lambda <- x$lambda.BIC
		pi <- x$pi.BIC
	}
	if(is.numeric(cluster.choice) == TRUE) {
		tmp <- names(x$BIC.all)
		clust.nums <- as.numeric(lapply(strsplit(tmp, "="), function(x) x[2]))
		index <- which(clust.nums == cluster.choice)
		if(length(index) == 0) {
			stop(paste(sQuote("cluster.choice"), sep = ""), " must be one of ",
				paste(dQuote("ICL"), sep = ""), " or ", paste(dQuote("BIC"), sep = ""), 
				" or \n one of the cluster sizes run for object ", paste(sQuote("x"), sep = ""))
		}
		probaPost <- x$probaPost[[index]]
		g <- clust.nums[[index]]
		labels <- x$labels[,index]
		lambda <- x$lambda[[index]]
		pi <- x$pi[[index]]
	}

	map <- apply(probaPost, 1, max)
	length(which(map > 0.9))/length(map)

	cat("*************************************************\n")
	cat("Selected number of clusters = ", g, sep = "")
	if(cluster.choice == "ICL") cat(" (ICL = ", x$ICL, ")\n", sep = "")
	if(cluster.choice == "BIC") cat(" (BIC = ", x$BIC, ")\n", sep = "")
	if(is.numeric(cluster.choice) == TRUE) cat("\n")
	cat("*************************************************\n")
	tab <- table(labels)
	names(tab) <- paste("Cluster", names(tab))
	cat("Cluster sizes:\n"); print(tab); cat("\n")
	cat("Number of observations with MAP > 0.90 (% of total):\n")
	cat(length(which(map > 0.9)), " (", round(length(which(map > 0.9))/length(map)*100,2),
		"%)\n\n", sep = "")
	cat("Number of observations with MAP > 0.90 per cluster (% of total per cluster):\n"); 

	tab2 <- matrix(NA, nrow = 2, ncol = g)
	colnames(tab2) <- names(tab); rownames(tab2) <- rep("", 2)
	for(i in 1:g) {
		if(sum(labels == i) > 1) {
			map.clust <- apply(probaPost[labels == i,], 1, max)
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 1) {
			map.clust <- max(probaPost[labels == i,])
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 0) {
			tab2[1,i] <- "---"
			tab2[2,i] <- "---"
		}
	}
	print(tab2, quote = FALSE); cat("\n")

	cat("Lambda:\n"); print(lambda); cat("\n")
	cat("Pi:\n"); print(pi); cat("\n")
}



