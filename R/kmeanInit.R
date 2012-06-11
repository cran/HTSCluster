kmeanInit <-
function(y, g, conds, lib.size, lib.type = "TC") {
if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0 | sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(g) != 1)
	stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
if(g < 0 | round(g) != g) 
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
n <- dim(y)[1];cols <- dim(y)[2];
y <- as.matrix(y, nrow = n, ncol = cols)
d <- length(unique(conds))
r <- as.vector(table(conds))
w <- rowSums(y)
if(lib.size == FALSE) {
s <- rep(1, cols)
}
if(lib.size == TRUE) {
if(lib.type == "TC") s <- colSums(y) / sum(y);
if(lib.type == "Q") s <- apply(y, 2, quantile, 0.75) / sum(apply(y, 2, quantile, 0.75));
if(lib.type == "MedRatio") {
tmp <- apply(y / (apply(y, 1, prod)^(1/cols)), 2, median, na.rm = TRUE)
s <- tmp/sum(tmp);
} 
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
pi.init <- as.vector(table(g.init$cluster)/n)
return(list(pi.init = pi.init, lambda.init = lambda.init))
}

