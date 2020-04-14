rm(list=ls())

library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()


## inverse cdf transform

densityq <- function(x) exp(dbeta(x, 4, 4, log=T) + cos(5*pi*x))
constantq <- integrate(densityq, 0, 1)$val

## cdf q
cdfq <- function(x) integrate(densityq, lower=0, upper=x)$val / constantq
## inverse cdf q
inversecdfq <- function(z){
    uniroot(function(x) cdfq(x) - z, interval=c(0, 1))$root
}
nsamples <- 1e4
u <- runif(nsamples, min=0, max=1)
x <- sapply(u, inversecdfq)


gscatter <- qplot(x=x, y=u, geom="blank") +
    geom_point(alpha=0.05)
gmargx <- qplot(x=x, geom="blank") +
    geom_histogram(aes(y=..density..), fill=graphsettings$colors[2]) +
    stat_function(fun= function(x) densityq(x) / constantq)
gmargy <- qplot(x=u, geom="blank") +
    geom_histogram(aes(y=..density..), fill=graphsettings$colors[1]) +
    stat_function(fun=dunif) + coord_flip() + scale_y_reverse()
empty <- ggplot()
g <- gridExtra::grid.arrange(empty, gmargx, gmargy, gscatter,
                             ncol=2, nrow=2,
                             widths=c(1, 4), heights=c(1, 4))
g
ggsave(filename="../icdf.pdf", plot=g)
