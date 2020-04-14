rm(list=ls())

library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()

## rejection sampling
nsamples <- 1e3
x <- mvtnorm::rmvnorm(nsamples,
                      mean=c(0, 0),
                      sigma=matrix(c(1, 0.8, 0.8, 1), nrow=2, ncol=2))
g <- qplot(x[, 1], x[, 2], geom="blank") +
    geom_point() +
    xlim(-4, 4) +
    ylim(-4, 4)
# constraint: x + y <= 1
rtrunc <- function(){
    accepted <- FALSE
    while (!accepted){
        x=mvtnorm::rmvnorm(1,
                           mean=c(0,0),
                           sigma=matrix(c(1, 0.8, 0.8, 1), nrow=2, ncol=2))
        if (x[1, 1] + x[1, 2] <= 1){
            return(x)
        }
    }
}
for (isample in 1:nsamples) x[isample,] <- rtrunc()
g <- qplot(x[, 1], x[, 2], geom="blank") +
    geom_point() +
    xlim(-4, 4) +
    ylim(-4, 4) +
    geom_abline(slope=-1, intercept=1)
g
## or...
## sample from a beta(3, 2) by sampling from uniforms

minimizer <- optimize(f=function(x) - dbeta(x, 3, 2, log=T),
                      interval=c(0, 1))$minimum
logM <- dbeta(minimizer, 3, 2, log=T)
#
y <- runif(nsamples)
us <- runif(nsamples)
accepts <- log(us) < (dbeta(y, 3, 2, log=TRUE) - 0 - logM)
x <- y[accepts]

pt.df <- data.frame(u=us, y=y, accept=accepts)
g <- ggplot(pt.df, aes(x=y, y=us  * exp(logM), colour=accept)) +
    geom_point() +
    scale_color_manual(values=graphsettings$colors) +
    stat_function(fun=function(x) dbeta(x, 3, 2), colour="black") +
    theme(legend.position="none")
g
ggsave(filename="../rejection.pdf", plot=g)
