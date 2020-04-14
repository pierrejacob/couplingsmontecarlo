rm(list=ls())

library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()

## importance sampling
rq <- function(n) rgamma(n, 2, 1)
logdq <- function(x) dgamma(x, 2, 1, log=TRUE)
logdp <- function(x) dgamma(x, 4, 2, log=TRUE)
nsamples <- 1e4
x <- rq(nsamples)
logw <- logdp(x) - logdq(x)
mlw <- max(logw)
mlw + log(mean(exp(logw - mlw)))

w <- exp(logw - mlw)
nw <- w / sum(w) # normalized weights

ghist1 <- qplot(x=x, geom="blank") +
    geom_histogram(aes(y=..density..),
                   alpha=1, fill=graphsettings$colors[1]) +
    stat_function(fun=function(x) exp(logdq(x)),
                  colour=graphsettings$colors[1]) +
    ylim(0,.6) +
    xlim(0, 12)
# geom_histogram(aes(y=..density.., weight=nw), alpha=0.75, fill=colors[2]) +
# stat_function(fun=function(x) exp(logdp(x)), colour=colors[2])

glogw <- qplot(x=x, xend=x, y=0, yend=logw, geom="blank") +
    geom_segment(alpha=1) +
    xlim(0, 12) +
    ylim(-7, 2)

ghist2 <- qplot(x=x, geom="blank") +
    geom_histogram(aes(y=..density.., weight=nw),
                   alpha=1, fill=graphsettings$colors[2]) +
    stat_function(fun=function(x) exp(logdp(x)),
                  colour=graphsettings$colors[2]) +
    ylim(0, .6) +
    xlim(0, 12)

g <- gridExtra::grid.arrange(ghist1, glogw, ghist2, ncol=1)
g
ggsave(filename="../importancesampling.pdf", plot=g)
