rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
set_theme_chapter3 <- function(){
    library(ggplot2)
    library(gridExtra)
    theme_set(theme_void())
    colors <- c(rgb(0.8,0.5,0.2), rgb(0.2, 0.6, 0.9))
    return(list(colors = colors))
}
graphsettings <- set_theme_chapter3()
library(ggridges)
library(reshape2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

## log target
# logtarget <- function(x){
#     evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = c(1.2,2.1), log = TRUE)
#     return(max(evals) + log(sum(exp(evals - max(evals)))))
# }

logtarget <- function(x) -abs(x)^3

## sd of proposal distribution for random walk MH
sd_proposal <- .5
## MH kernel
MH_kernel <- function(chain_state){
    proposal_value <- rnorm(1, chain_state, sd_proposal)
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    if (log(runif(1)) < (proposal_pdf - current_pdf)){
        return(proposal_value)
    } else {
        return(chain_state)
    }
}
##

## target
curve(sapply(x, function(xp) exp(logtarget(xp))), from = -4, to = 4)
## possible drift function
## curve(sapply(x, function(xp) exp(-0.5*logtarget(xp))), col = "red", from = -5, to = 8)

## compute rejection probability for each point by numerical integration
rejectproba <- function(x0){
    logtargetx <- logtarget(x0)
    integrand_ <- function(y){
        sapply(y, function(v) pmin(1, exp(logtarget(v) - logtargetx)) *
                   dnorm(v, mean = x0, sd = sd_proposal))
    }
    return(1-integrate(integrand_, lower = x0 - 5*sd_proposal, upper = x0 + 5*sd_proposal, subdivisions = 1e3)$val)
}
curve(sapply(x, rejectproba), from = -15, to = 8, n = 1e3, ylim = c(0,1))

continuouspart <- function(x0){
    rejecx0 <- rejectproba(x0)
    logtargetx <- logtarget(x0)
    f <- function(y) (1-rejecx0) * pmin(1, exp(logtarget(y) - logtargetx)) * dnorm(y, mean = x0, sd = sd_proposal)
    return(f)
}

xseq <- seq(from = -6, to = 6, length.out = 500)
continuouspart.df <- data.frame()
plot(NULL, xlim = c(-6,6), ylim = c(0,.8), xlab = "x", ylab = "density")
for (x0 in seq(from = -4, to = 4, length.out = 50)){
    continuouspart_x0 <- continuouspart(x0)
    curve(sapply(x, continuouspart_x0), add = TRUE, n = 1e3)
    # abline(v = x0, lty = 3)
    continuouspart.df <- rbind(continuouspart.df, data.frame(x = xseq, y = sapply(xseq, continuouspart_x0),
                                                             x0 = x0))
}

##

##
g <- ggplot(continuouspart.df, aes(x = x, y = y, group = x0)) + geom_line(colour = graphsettings$colors[1])
# g + stat_function(fun = function(z) exp(logtarget(z)))
g
# ggsave(filename="../minorization.pdf", plot = g, width = 10, height = 5)

plot(NULL, xlim = c(-6,6), ylim = c(-150,0), xlab = "x", ylab = "log-density")
for (x0 in seq(from = -4, to = 4, length.out = 25)){
    continuouspart_x0 <- continuouspart(x0)
    curve(log(sapply(x, continuouspart_x0)), add = TRUE, n = 1e3)
    # abline(v = x0, lty = 3)
}

##
# V <- function(x) exp(-0.5*logtarget(x))
V <- function(x) (1+x^2)
## expected value of V given x0
expectednextV <- function(x0){
    rejecx0 <- rejectproba(x0)
    logtargetx <- logtarget(x0)
    f <- function(y) (1-rejecx0) * pmin(1, exp(logtarget(y) - logtargetx)) *
        dnorm(y, mean = x0, sd = sd_proposal) *
        V(y)
    result <- V(x0) * rejecx0 + integrate(f, lower = x0 - 5*sd_proposal, upper = x0 + 5*sd_proposal)$val
    return(result)
}
curve(sapply(x, expectednextV), from = -4, to = 4)
curve(0.8 * sapply(x, V), add = T, col = 'red')

curve(sapply(x, expectednextV)/sapply(x, V), from = -4, to = 4)

curve(log(sapply(x, expectednextV)) - log(sapply(x, V)), from = -10, to = 10)
abline(h = 0)

