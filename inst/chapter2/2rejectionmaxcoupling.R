rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter2()

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

p <- function(x) dbeta(x, 5, 5)
q <- function(x) dbeta(x, 2, 1.5)

rmc <- getmaxcoupling(rp=function(n) rbeta(n, 5, 5),
                     dp=function(x) dbeta(x, 5, 5, log=TRUE),
                     rq=function(n) rbeta(n, 2, 1.5),
                     dq=function(x) dbeta(x, 2, 1.5, log=TRUE))

res <- sapply(1:1e4, function(i) rmc()$xy)

gscatter <- qplot(x=res[1,], y=res[2,],
                 geom="blank") +
            geom_point(alpha=0.5)
gmargx <- qplot(x=res[1,], geom="blank") +
            geom_histogram(aes(y=..density..),
                             fill=graphsettings$colors[1])  +
            stat_function(fun=p)
gmargy <- qplot(x=res[2,], geom="blank") +
            geom_histogram(aes(y=..density..), fill=graphsettings$colors[2]) +
            stat_function(fun=q) +
            coord_flip()
empty <- ggplot()

g <- gridExtra::grid.arrange(gmargx, empty, gscatter, gmargy,
                            ncol=2, nrow=2,
                            widths=c(4, 1), heights=c(1, 4))

ggsave(filename="../tvdistance2.pdf", plot = g)
