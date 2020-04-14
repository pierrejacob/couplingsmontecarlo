## snippets of code
## compute average of weights on log-scale
N <- 1e4
x <- rnorm(N)
logf <- function(x) dcauchy(x, log=TRUE)
logg <- function(x) dnorm(x, log=TRUE)
logw <- logf(x) - logg(x)
mlw <- max(logw)
mlw + log(mean(exp(logw - mlw)))

## with some probability do this, otherwise do that
if (log(runif(1)) < (logf(x[1]) - logg(x[1]))){
    Z <- -Inf
} else {
    Z <- x
}

## sample a bivariate Markov chain
T <- 50
x <- matrix(0, nrow=T, ncol=2)
xcurrent <- c(1,-1)
for (t in 1:T){
    increment <- rnorm(1)
    xcurrent <- 0.9 * xcurrent + increment
    x[t,] <- xcurrent
}
