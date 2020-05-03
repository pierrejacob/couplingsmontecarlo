## rejection sampler to obtain pairs (X,Y)
## such that X ~ p and Y ~ q
## and {X=Y} occurs with maximal probability
## rp: function to sample from p
## dp: function to evaluate the log-density of p
## rq: function to sample from q
## dq: function to evaluate the log-density of q
#'@export
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
