summary.HTSClusterWrapper <-
function (object, ...) 
{
	x <- object
    	if (class(x) != "HTSClusterWrapper") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSClusterWrapper"), sep = ""), sep = "")
              }
        
	cat("*************************************************\n")
	cat("Selected number of clusters = ", ncol(x$select.results$lambda), "\n", sep = "")
	cat(" (ICL = ", x$select.results$ICL, ")\n", sep = "")
	cat(" (BIC = ", x$select.results$BIC, ")\n", sep = "")
	cat("*************************************************\n")
        cat(" ICL values for all fitted models: \n")
        print(x$ICL.all)
	cat("*************************************************\n")

      }

