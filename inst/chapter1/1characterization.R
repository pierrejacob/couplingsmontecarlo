rm(list=ls())

library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()

## two distributions in [0,1]
## density p
densityp <- function(x){
    0.25 * dbeta(x, 40,10) +
        0.5 * dbeta(x, 20,20) +
        0.25 * dbeta(x, 10, 45)
}
## density q
densityq <- function(x) exp(dbeta(x, 4, 4, log=T) + cos(5*pi*x))
constantq <- integrate(densityq, 0, 1)$val

## cdf p
cdfp <- function(x) integrate(densityp, lower=0, upper=x)$val


## inverse cdf p
inversecdfp <- function(z){
    uniroot(function(x) cdfp(x) - z, interval=c(0, 1))$root
}

xseq <- seq(from = 0, to = 1, length.out = 1e3)

gdensity <- qplot(xseq, sapply(xseq, densityp), geom = "line")
gcdf <- qplot(xseq, sapply(xseq, cdfp), geom = "line")
ggrad <- qplot(xseq[40:960], sapply(xseq[40:960], function(z) numDeriv::grad(function(t) log(densityp(t)), z)), geom = "line", xlim = c(0,1))
ghist <- qplot(x = sapply(runif(1e4), inversecdfp), geom = "blank") + geom_histogram(aes(y=..density..), fill = graphsettings$colors[2])


g <- gridExtra::grid.arrange(gdensity, gcdf, ggrad, ghist,
                             ncol=2, nrow=2)
g
ggsave(filename = "../characterizations.pdf", plot = g)
