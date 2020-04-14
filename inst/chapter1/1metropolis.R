rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()


# target log-density evaluation
bananatarget <- function(x) -(1-x[1])^2 - 10*((x[2]-x[1]^2)^2)
# gradient of target log-density
bananagradtarget <- function(x) c(-2*(x[1]-1) + 40*x[1]*(x[2]-x[1]^2),
                                  -20 * (x[2]-x[1]^2))

logdtarget <- bananatarget
logdgradient <- bananagradtarget

nmcmc <- 1000
dimension <- 2
## initialization
state <- list(x=c(4, -4))
state$logdtarget <- logdtarget(state$x)
## store entire history
chain <- matrix(nrow=dimension, ncol=nmcmc)
proposal_history <- matrix(nrow=dimension, ncol=nmcmc)
## number of accepted proposals
naccepts <- 0
for (imcmc in 1:nmcmc){
    ## propose new point
    proposal <- state$x + rnorm(dimension, mean=0, sd=1)
    proposal_history[,imcmc] <- proposal
    ## compute target density at proposed point
    proposal_logdtarget <- logdtarget(proposal)
    ## decide to accept or not
    accept <- log(runif(1)) < proposal_logdtarget - state$logdtarget
    if (accept){
        ## if accept
        naccepts <- naccepts + 1
        state <- list(x=proposal, logdtarget=proposal_logdtarget)
    }
    ## else, no modification to state
    ## store current state
    chain[, imcmc] <- state$x
}
cat("accept rate", 100 * naccepts / nmcmc, "%\n")

matplot(t(chain), type = 'l')

plot(chain[1,], chain[2,])
## plot path on 2d space

chaindf <- data.frame(x = chain[1,1:(nmcmc-1)],
           xend = chain[1,2:(nmcmc)],
           y = chain[2,1:(nmcmc-1)],
           yend = chain[2,2:(nmcmc)])
chaindf$iteration <- 1:(nmcmc-1)
head(chaindf)
ggplot(chaindf, aes(x=x, xend=xend, y=y, yend=yend)) + geom_path(arrow = arrow(length = unit(.25, 'inches'))) + geom_point()


# ## 3d trace plot
# x <- c(-2,3)
# yrange <- c(-1,6)
# z <- matrix(rep(0, 4), ncol = 2)
# op <- par(bg = "white")
# # png(filename = "3dtrace.png", height = 600, width = 600, pointsize = 12)
# persp_res <- persp(x, yrange, z, theta = 60, phi = 30, d = 2, r = 3, expand = 1, col = 'white',
#                    zlim = c(0,1000), border = NA, box = FALSE, axes = 0)
#
# lines(trans3d(chaindf$x, chaindf$y, chaindf$iteration, pmat = persp_res), col = rgb(0,0,0))

## MALA algorithm

h <- 0.005
nmcmc <- 1000
dimension <- 2
## initialization
state <- list(x=c(4, -4))
state$logdtarget <- logdtarget(state$x)
## store entire history
malachain <- matrix(nrow=dimension, ncol=nmcmc)
malaproposal_history <- matrix(nrow=dimension, ncol=nmcmc)
## number of accepted proposals
naccepts <- 0
for (imcmc in 1:nmcmc){
    ## propose new point
    proposal <- state$x + (h/2) * logdgradient(state$x) + rnorm(dimension, mean=0, sd=sqrt(h))
    malaproposal_history[,imcmc] <- proposal
    ## compute target density at proposed point
    proposal_logdtarget <- logdtarget(proposal)
    ## decide to accept or not
    logmhaccept <- (proposal_logdtarget - state$logdtarget)
    logmhaccept <- logmhaccept + sum(dnorm(state$x,  proposal + (h/2) * logdgradient(proposal), sd = sqrt(h), log=TRUE))
    logmhaccept <- logmhaccept - sum(dnorm(proposal, state$x + (h/2) * logdgradient(state$x),   sd = sqrt(h), log=TRUE))
    accept <- log(runif(1)) < logmhaccept
    if (accept){
        ## if accept
        naccepts <- naccepts + 1
        state <- list(x=proposal, logdtarget=proposal_logdtarget)
    }
    ## else, no modification to state
    ## store current state
    malachain[, imcmc] <- state$x
}
cat("accept rate", 100 * naccepts / nmcmc, "%\n")

matplot(t(malachain), type = 'l')

plot(malachain[1,], malachain[2,])
## plot path on 2d space



# ## 3d trace plot
# x <- c(-2,3)
# yrange <- c(-1,6)
# z <- matrix(rep(0, 4), ncol = 2)
# op <- par(bg = "white")
# # png(filename = "3dtrace.png", height = 600, width = 600, pointsize = 12)
# persp_res <- persp(x, yrange, z, theta = 60, phi = 30, d = 2, r = 3, expand = 1, col = 'white',
#                    zlim = c(0,1000), border = NA, box = FALSE, axes = 0)
#
# lines(trans3d(chaindf$x, chaindf$y, chaindf$iteration, pmat = persp_res), col = rgb(0,0,0))







