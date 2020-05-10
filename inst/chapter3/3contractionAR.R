rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)

graphsettings <- set_theme_chapter3()
library(ggridges)
library(reshape2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

rho <- 0.95
ar_function <- function(x, u) rho*x + u

## forward chain
x <- rnorm(1, 0, 10)
y <- rnorm(1, 0, 10)
nmcmc <- 200
forwardxchain <- rep(0, nmcmc)
forwardychain <- rep(0, nmcmc)
us <- rep(0, nmcmc)
for (imcmc in 1:nmcmc){
    u <- rnorm(1, 0, 1)
    x <- ar_function(x, u)
    y <- ar_function(y, u)
    forwardxchain[imcmc] <- x
    forwardychain[imcmc] <- y
    us[imcmc] <- u
}

matplot(cbind(forwardxchain, forwardychain), type = 'l')

## backward process
x0 <- rnorm(1,0,10)
y0 <- rnorm(1,0,10)
backwardxchain <- rep(0, nmcmc)
backwardychain <- rep(0, nmcmc)
xalt <- rep(0, 1+nmcmc)
xalt[1] <- x0
yalt <- rep(0, 1+nmcmc)
for (imcmc in 1:nmcmc){
    xcurrent <- x0
    ycurrent <- y0
    for (j in 1:imcmc){
        xcurrent <- ar_function(xcurrent, us[imcmc-j+1])
        ycurrent <- ar_function(ycurrent, us[imcmc-j+1])
    }
    backwardxchain[imcmc] <- xcurrent
    backwardychain[imcmc] <- ycurrent
    xalt[imcmc+1] <- xalt[imcmc] + (rho^(imcmc-1))*(us[imcmc] + (rho - 1)*xalt[1])
}

matplot(cbind(backwardxchain, backwardychain), type = 'l')
## matplot(cbind(xchain, xalt[2:(nmcmc+1)]), type = 'l')
#

gforward <- ggplot(rbind(data.frame(time = 1:nmcmc, chain = forwardxchain, i = 1),
                         data.frame(time = 1:nmcmc, chain = forwardychain, i = 2)),
                         aes(x = time, y = chain, group = factor(i), colour = factor(i))) + geom_line()
gforward <- gforward + theme(legend.position = "none") + scale_color_manual(values = graphsettings$colors)

gbackward <- ggplot(rbind(data.frame(time = 1:nmcmc, chain = backwardxchain, i = 1),
                          data.frame(time = 1:nmcmc, chain = backwardychain, i = 2)),
                   aes(x = time, y = chain, group = factor(i), colour = factor(i))) + geom_line()
gbackward <- gbackward + theme(legend.position = "none") + scale_color_manual(values = graphsettings$colors)
g <- gridExtra::grid.arrange(gforward, gbackward, nrow = 2)
g
ggsave(filename="../forwardbackward.pdf", plot = g, width = 10, height = 7)
