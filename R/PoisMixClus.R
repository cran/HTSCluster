PoisMixClus <-
function(y, gmin, gmax, lib.size = TRUE, lib.type = "TC", conds, 
init.type = "small-em", alg.type = "EM", cutoff = 10e-6, iter = 1000, mean.filter = FALSE,
verbose = FALSE) {

if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
	stop(paste(sQuote("y"), "must be a matrix"))
if(min(y) < 0 | sum(round(y)) != sum(y)) 
	stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
if(min(rowSums(y)) == 0 & mean.filter == FALSE)
	stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
if(length(gmin) != 1)
	stop(paste(sQuote("gmin"), "(the minimum number of clusters) must be a nonnegative integer"))
if(gmin < 0 | round(gmin) != gmin) 
	stop(paste(sQuote("gmin"), "(the minimum number of clusters) must be a nonnegative integer"))
if(length(gmax) != 1)
	stop(paste(sQuote("gmax"), "(the maximum number of clusters) must be a nonnegative integer"))
if(gmax < 0 | round(gmax) != gmax) 
	stop(paste(sQuote("gmax"), "(the maximum number of clusters) must be a nonnegative integer"))
if(gmin > gmax)
	stop(paste(sQuote("gmin"), "must be less than or equal to", sQuote("gmax")))
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
if(length(init.type) > 1)
	stop(paste(sQuote("init.type"), "must be of length 1"))
if(init.type != "small-em" & init.type != "kmeans") 
	stop(paste(sQuote("init.type"), "must be one of", dQuote("small-em"), "or", dQuote("kmeans")))
if(alg.type != "EM" & alg.type != "CEM")
	stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
if(length(alg.type) > 1)
	stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
if(is.logical(verbose) == FALSE)
	stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))

## Grouping columns of y in order of condition (all replicates put together)
o.ycols <- order(conds)
y <- y[,o.ycols]
conds <- conds[o.ycols]
conds.names <- unique(conds)

d <- length(unique(conds))
r <- as.vector(table(conds))
diff <- 100 ## Convergence criterion
if(length(rownames(y)) == 0) rn <- 1:nrow(y);
if(length(rownames(y)) > 0) rn <- rownames(y);
y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
rownames(y) <- rn;
## Remove genes with very low overall counts
ind.remove <- NA
if(mean.filter != FALSE & is.numeric(mean.filter) == TRUE) {
	ind.remove <- which(rowMeans(y) < mean.filter)
	if(length(ind.remove) > 0) y <- y[-ind.remove,];
}
n <- dim(y)[1];cols <- dim(y)[2]

w <- rowSums(y)
s <- rep(NA, cols)
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

BIC <- ICL <- as.data.frame(matrix(NA, nrow = 1, ncol = gmax-gmin+1))
colnames(BIC) <- colnames(ICL) <- paste("g=",gmin:gmax,sep="")
labels <- as.data.frame(matrix(NA, nrow = n, ncol = gmax-gmin+ 1))
colnames(labels) <- paste("g=",gmin:gmax,sep="")
probaPost <- lambda.final <- pi.final <- vector("list", gmax - gmin + 1)
K.index <- 1; g.values <- gmin:gmax

## Choose value of g
for(K in gmin:gmax) {

index <- 0;go <- 1;

## Inital values
## init.type: "kmeans", "small-em"
init.args <- list(y, K, conds, lib.size)
if(init.type == "kmeans") init.alg <- "kmeanInit";
if(init.type == "small-em") {
init.alg <- "emInit"
init.args <- list(y = y, g = K, conds = conds, lib.size = lib.size, 
lib.type = lib.type, alg.type = alg.type, starts = 5, verbose = verbose)
}
param.init <- do.call(init.alg, init.args)

pi <- pi.old <- param.init$pi.init
lambda <- lambda.old <- param.init$lambda.init
mean.calc <- mean.old <- PoisMixMean(y = y, g = K, conds = conds, 
	s = s, lambda = lambda)

while(go == 1) {

############
## E-step ##
############
t <- probaPost(y, K, conds, pi, s, lambda)

############
## C-step ##
############
if(alg.type == "CEM") {
## If two values of t_{ik} are map, 
## arbitrarily choose the first
partition <- unlist(apply(t, 1, 
function(x) which(x == max(x, na.rm = TRUE))[1]))
partition.mat <- matrix(0, nrow = n, ncol = K)
for(i in 1:n) partition.mat[i,partition[i]] <- 1;
}

############
## M-step ##
############
if(alg.type == "CEM") {
for(k in 1:K) {
pi[k] <- length(which(partition == k))/n
}
denom <- colSums(partition.mat * w)
for(j in 1:d) {
denom.bis <- denom * s.dot[j]
num <- colSums(partition.mat * 
rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])))
lambda[j,] <- num / denom.bis
}
}

if(alg.type == "EM") {
pi <- colSums(t)/n
denom <- colSums(t * w)
for(j in 1:d) {
denom.bis <- denom * s.dot[j]
num <- colSums(t * 
matrix(rep(rowSums(as.matrix(y[,which(conds == (unique(conds))[j])])),K), 
ncol = K))
lambda[j,] <- num / denom.bis
}
}

#################
## Convergence ##
#################
mean.calc <- PoisMixMean(y, g = K, conds, s, lambda)
diff <- abs(logLikePoisMixDiff(y, mean.calc, pi, mean.old, pi.old))
lambda.old <- lambda; pi.old <- pi; mean.old <- mean.calc;

index <- index + 1
if(verbose == TRUE) print(paste("Log like diff:", diff))
if(diff < cutoff) go <- 0;
if(iter != FALSE & iter == index) go <- 0;
}

#####################################
## Final estimates of lambda and p ##
#####################################
names(pi) <- paste("Cluster", 1:K)
colnames(lambda) <- paste("Cluster", 1:K)
rownames(lambda) <- conds.names
lambda.final[[K.index]] <- lambda
pi.final[[K.index]] <- pi

## Check to make sure one of the components is not degenerate
if(min(pi) == 0 | is.nan(sum(lambda)) == TRUE) {
probaPost[[K.index]] <- NA
labels[[K.index]] <- NA
BIC[,K.index] <- NA
ICL[,K.index] <- NA
}
if(min(pi) > 0 | is.nan(sum(lambda)) == FALSE) {

mean.calc <- PoisMixMean(y, K, conds, s, lambda)
LL.tmp <- logLikePoisMix(y, mean.calc, pi)
LL <- LL.tmp$ll

######################
## Determine labels ##
######################
t <- probaPost(y, K, conds, pi, s, lambda)
## If two clusters have exactly identical map estimators,
## arbitrarily choose the first one
map <- unlist(apply(t, 1, function(x) which(x == max(x, 
na.rm = TRUE))[1]))
z <- matrix(0, nrow = n, ncol = K)
for(i in 1:n) z[i,map[i]] <- 1;
probaPost[[K.index]] <- t
labels[[K.index]] <- map

##############################
## Calculate BIC, ICL, SICL ##
##############################
np <- (K-1) + n + (d-1)*K # pi + w + lambda
BIC[,K.index] <- -LL + (np/2) * log(n)
entropy <- -2*sum(z*log(t), na.rm = TRUE)
ICL[,K.index] <- BIC[,K.index] + entropy
}
K.index <- K.index + 1
}

## CHOSEN MODELS BY BIC AND ICL
## RETURN FITTED VALUES AND RESIDUALS
BIC.choice <- ICL.choice <- g.BICchoice <- g.ICLchoice <- 
lab.BIC <- lab.ICL <- lambda.BIC <- lambda.ICL <- pi.BIC <- pi.ICL <- 
probaPost.BIC <- probaPost.ICL <- NA

if(sum(is.na(BIC) == TRUE) < length(BIC)) {
if(min(BIC, na.rm = TRUE) < Inf) {
BIC.choice <- min(BIC, na.rm = TRUE)
ICL.choice <- min(ICL, na.rm = TRUE)
g.BICchoice <- c(gmin:gmax)[which(BIC == BIC.choice)]
g.ICLchoice <- c(gmin:gmax)[which(ICL == ICL.choice)]
lab.BIC <- labels[[which(BIC == BIC.choice)]]
lab.ICL <- labels[[which(ICL == ICL.choice)]]
probaPost.BIC <- probaPost[[which(BIC == BIC.choice)]]
probaPost.ICL <- probaPost[[which(ICL == ICL.choice)]]
lambda.BIC <- lambda.final[[which(BIC == BIC.choice)]]
lambda.ICL <- lambda.final[[which(ICL == ICL.choice)]]
pi.BIC <- pi.final[[which(BIC == BIC.choice)]]
pi.ICL <- pi.final[[which(ICL == ICL.choice)]]

}
}

results <- list(lambda = lambda.final, pi = pi.final, labels = labels, 
probaPost = probaPost, BIC.all = -BIC, ICL.all = -ICL, 
alg.type = alg.type, BIC = -BIC.choice, ICL = -ICL.choice, 
g.BIC = g.BICchoice, g.ICL = g.ICLchoice, labels.BIC = lab.BIC,
labels.ICL = lab.ICL, lambda.BIC = lambda.BIC, pi.BIC = pi.BIC,
lambda.ICL = lambda.ICL, pi.ICL = pi.ICL, probaPost.BIC = probaPost.BIC,
probaPost.ICL = probaPost.ICL, lib.size = lib.size, lib.type = lib.type, s = s,
y = y, conds = conds, ind.remove = ind.remove, mean.filter = mean.filter)

class(results) <- "HTSCluster"
return(results)
}

