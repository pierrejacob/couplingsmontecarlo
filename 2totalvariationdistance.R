rm(list=ls())

library(ggplot2)
theme_set(theme_void())

colors <- c(rgb(1, 0.1, 0.3), rgb(0.3, 0.1, 1))

## visualization of TV bounds
p <- function(x) dbeta(x, 5, 5)
q <- function(x) dbeta(x, 2, 1.5)

xseq <- seq(from=0, to=1, length.out=1e3)
pseq <- p(xseq)
qseq <- q(xseq)

g1 <- ggplot(data=data.frame(x=xseq, y=pseq), aes(x=x, y=y)) +
        geom_polygon(fill=colors[1], alpha=0.75) +
        geom_polygon(aes(y=qseq),
                       fill=colors[2], alpha=0.75) +
        geom_polygon(aes(y=pmin(pseq,qseq)),
                       col='black', alpha=0.5, linetype=2) +
                            xlim(0,1)
g1

g2 <- ggplot(data=data.frame(x=xseq, y=pseq-qseq), aes(x=x, y=y)) +
        geom_polygon(fill=rgb(0.9,0.8,0.2), alpha=1) +
        stat_function(fun=function(x){ 1*(p(x)>q(x)) - 1*(p(x)<q(x)) }) +
        xlim(0,1)
g2

