## This script contains general definitions used in all other scripts

## load Rcpp
library(Rcpp)

## graphics
library(ggplot2)
colors <- c(rgb(1, 0.1, 0.3), rgb(0.3, 0.1, 1))
theme_set(theme_void())

## rejection sampler to obtain pairs (X,Y)
## such that X ~ p and Y ~ q
## and {X=Y} occurs with maximal probability
## rp: function to sample from p
## dp: function to evaluate the log-density of p
## rq: function to sample from q
## dq: function to evaluate the log-density of q
getmaxcoupling <- function(rp, dp, rq, dq){
	function() {
		x <- rp(1)
		if (dp(x) + log(runif(1)) < dq(x)) {
			return(list(xy=c(x, x), identical=TRUE))
		}
		else {
			reject <- TRUE
			y <- NA
			while (reject) {
				y <- rq(1)
				reject <- (dq(y) + log(runif(1)) < dp(y))
			}
			return(list(xy=c(x, y), identical=FALSE))
		}
	}
}
