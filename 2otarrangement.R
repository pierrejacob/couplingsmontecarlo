rm(list = ls())
source("commongrounds.R")
library(doParallel)
library(doRNG)
registerDoParallel(cores = 10)

## two distributions in [0,1]
## density p
densityp <- function(x){ 0.25 * dbeta(x, 40,10) + 0.5 * dbeta(x, 20,20) + 0.25 * dbeta(x, 10, 45) }
## density q
densityq <- function(x){exp(dbeta(x, 4, 4, log = T) + cos(5*pi*x)) }
constantq <- integrate(densityq, 0, 1)$val

## cdf p
cdfp <- function(x) integrate(densityp, lower = 0, upper = x)$val
## cdf q
cdfq <- function(x) integrate(densityq, lower = 0, upper = x)$val / constantq

## inverse cdf p
inversecdfp <- function(z) uniroot(function(x) cdfp(x) - z, interval = c(0, 1))$root
## inverse cdf q
inversecdfq <- function(z) uniroot(function(x) cdfq(x) - z, interval = c(0, 1))$root

nsamples <- 1e4
unifs <- runif(nsamples)
xys <- foreach(irep = 1:nsamples) %dorng% c(inversecdfp(unifs[irep]), inversecdfq(unifs[irep]))
xys <- sapply(xys, function(x) x)
xs <- xys[1,]
ys <- xys[2,]

gsegments <- qplot(x = 0, xend = 1, y = xs, yend = ys, geom = "blank") +geom_segment(alpha = 0.05)
gmargleft <- qplot(x = xs, geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[1]) + stat_function(fun = function(x) densityp(x)) + coord_flip() + scale_y_reverse()
gmargright <- qplot(x = ys, geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[2]) + stat_function(fun = function(x) densityq(x)/constantq)  + coord_flip()
g <- gridExtra::grid.arrange(gmargleft, gsegments, gmargright, ncol=3, nrow=1, widths=c(1,4,1), heights=4)

# ggsave(filename = "otarrangementwithmarginals.pdf", plot = g)
