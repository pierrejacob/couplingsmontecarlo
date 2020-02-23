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


## two cdf
xseq = seq(from = 1e-10, to = 1-1e-10, length.out = 1000)
cdfp_seq <- sapply(xseq, cdfp)
cdfq_seq <- sapply(xseq, cdfq)
gcdf <- qplot(x = xseq, geom = "blank") + geom_line(aes(y = cdfp_seq), color = colors[1]) +
  geom_line(aes(y = cdfq_seq), color = colors[2])
gcdf

## two inverse cdf
invcdfp_seq <- sapply(xseq, inversecdfp)
invcdfq_seq <- sapply(xseq, inversecdfq)
ginvcdf <- qplot(x = xseq, geom = "blank") + geom_line(aes(y = invcdfp_seq), color = colors[1]) +
  geom_line(aes(y = invcdfq_seq), color = colors[2])
ginvcdf

## optimal transport map
otmap <- function(x) inversecdfq(cdfp(x))
otmap_seq <- sapply(xseq, otmap)
gotmap <- qplot(x = xseq, geom = "blank") + geom_line(aes(y = otmap_seq))
gotmap

gridExtra::grid.arrange(gcdf, ginvcdf, gotmap, nrow = 1)

# ggsave(filename = "cdfinverseotmap.pdf", plot = gridExtra::grid.arrange(gcdf, ginvcdf, gotmap, nrow = 1),
#        width = 15, height = 5)
