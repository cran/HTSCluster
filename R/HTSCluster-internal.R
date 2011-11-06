.d0Func <-
function(x, mean) {
return(ifelse(x == 0, mean, x * log(x/mean) + mean - x))
}

