plot.HTSClusterWrapper <-
function (x, file.name = FALSE, 
graphs = c("ICL", "BIC"), ...) 
{	
    	if (class(x) != "HTSClusterWrapper") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSClusterWrapper"), sep = ""), sep = "")
    	}
        
	if(file.name != FALSE) pdf(paste(file.name));
        
        if("ICL" %in% graphs & "BIC" %in% graphs) {
          par(mfrow = c(1,2), mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
          lines(gpl, x$ICL.all, lwd=2)
          points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
          plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", ylab = "BIC",
               main="BIC", pch=19)
          lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
        }

        if("ICL" %in% graphs & !"BIC" %in% graphs) {
          par(mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, x$ICL.all, xlab = "Number of clusters", ylab = "ICL", main="ICL", pch=19)
          lines(gpl, x$ICL.all, lwd=2)
          points(gpl[which(x$ICL.all == max(x$ICL.all))], max(x$ICL.all), col="red", pch="X", font=2, cex=2)
        }

        if(!"ICL" %in% graphs & "BIC" %in% graphs) {
          par(mar = c(4,4,2,2))
          gpl <- unlist(lapply(strsplit(names(x$ICL.all), split="=", fixed=TRUE), function(xx) xx[2]))
          plot(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), xlab = "Number of clusters", ylab = "BIC",
               main="BIC", pch=19)
          lines(gpl, unlist(lapply(x$all.results, function(xx) xx$BIC)), lwd=2)          
        }

	if(file.name != FALSE)  dev.off();
      }

