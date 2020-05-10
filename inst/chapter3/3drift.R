## illustration of drift function following Roberts and Rosenthal 2004
## Example 4 of that paper: Exponential target, random walk with uniform proposal on a small interval
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

logtarget <- function(x) dexp(x, log = TRUE)

## length of proposal interval for random walk MH
l_proposal <- 1
## MH kernel
MH_kernel <- function(chain_state){
    proposal_value <- runif(1, min = chain_state-l_proposal, max = chain_state+l_proposal)
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    if (log(runif(1)) < (proposal_pdf - current_pdf)){
        return(proposal_value)
    } else {
        return(chain_state)
    }
}

## compute rejection probability for each point by numerical integration
## this is very bruteforce; integrals could be computed exactly here
rejectproba <- function(x0){
    logtargetx <- logtarget(x0)
    integrand_ <- function(y){
        sapply(y, function(v) pmin(1, exp(logtarget(v) - logtargetx)) *
                   dunif(v, min = x0 - l_proposal, max = x0 + l_proposal))
    }
    return(1-integrate(integrand_, lower = x0 - l_proposal, upper = x0 + l_proposal, subdivisions = 1e3)$val)
}

continuouspart <- function(x0){
    rejecx0 <- rejectproba(x0)
    logtargetx <- logtarget(x0)
    f <- function(y) (1-rejecx0) * pmin(1, exp(logtarget(y) - logtargetx)) * dunif(y, min = x0 - l_proposal, max = x0 + l_proposal)
    return(f)
}

curve(sapply(x, rejectproba), from = 0, to = 2, n = 1e3, ylim = c(0,1))
plot(NULL, xlim = c(0,6), ylim = c(0,10), xlab = "x", ylab = "density")
for (x0 in seq(from = 0, to = 6, length.out = 50)){
    continuouspart_x0 <- continuouspart(x0)
    curve(sapply(x, continuouspart_x0), add = TRUE, n = 1e2)
    # abline(v = x0, lty = 3)
}


V <- function(x) exp(0.1*x)
## expected value of V given x0
expectednextV <- function(x0){
    rejecx0 <- rejectproba(x0)
    logtargetx <- logtarget(x0)
    f <- function(y) (1-rejecx0) * pmin(1, exp(logtarget(y) - logtargetx)) *
        dunif(y, min = x0 - l_proposal, max = x0 + l_proposal) *
        V(y)
    result <- V(x0) * rejecx0 + integrate(f, lower = x0 - l_proposal, upper = x0 + l_proposal)$val
    return(result)
}
curve(sapply(x, expectednextV), from = 0, to = 4, n = 1e3)
curve(sapply(x, V), add = T, col = 'red')
abline(v = l_proposal)

curve(log(sapply(x, expectednextV)), from = 0, to = 4, n = 1e3)
curve(log(sapply(x, V)), add = T, col = 'red')
abline(v = l_proposal)

curve(sapply(x, expectednextV)/sapply(x, V), from = 0, to = 4, ylim = c(0.7,1), n = 1e3)
abline(v = l_proposal)

##
independent_chains_from_x <- function(l_proposal, x, niteration, nrep = 100){
    chains <- matrix(NA, nrow = niteration+1, ncol = nrep)
    current_chains <- rep(x, nrep)
    current_pdf <- logtarget(current_chains)
    chains[1,] <- current_chains
    for (iteration in 1:niteration){
        proposal_values <- runif(nrep, min = current_chains-l_proposal, max = current_chains+l_proposal)
        proposal_pdf <- logtarget(proposal_values)
        accepts <- (log(runif(nrep)) < (proposal_pdf - current_pdf))
        current_chains[accepts] <- proposal_values[accepts]
        current_pdf[accepts] <- proposal_pdf[accepts]
        chains[iteration+1,] <- current_chains
    }
    return(chains)
}

chains <- independent_chains_from_x(l_proposal, x = 8, niteration = 5e3, nrep = 1000)
## find hitting times
hitting_times <- apply(chains, 2, function(v) which(v<l_proposal)[1])
hist(hitting_times, nclass = 50)
plot(NULL, xlim = c(0,200), ylim = c(0, max(chains)))
for (ichain in 1:250){
    lines(x = 1:hitting_times[ichain], y = chains[1:hitting_times[ichain],ichain], col = rgb(0,0,0,alpha=0.1))
}
points(x = hitting_times[1:250], y = sapply(1:250, function(i) chains[hitting_times[i], i]))
abline(h = l_proposal)
##
chaindf <- data.frame()
for (ichain in 1:ncol(chains)){
    chaindf <- rbind(chaindf,
                     data.frame(x = 1:hitting_times[ichain],
                                y = chains[1:hitting_times[ichain], ichain],
                                ichain = ichain))
}
head(chaindf)


chaindf <- chaindf %>% group_by(ichain) %>% mutate(noise = cos(2*pi*.01 + ichain) + rnorm(n())) %>%
    mutate(noise2 = (noise - min(noise))/(max(noise)-min(noise))) %>% ungroup()
g <- ggplot(chaindf, aes(x = x, y = y, group = ichain, alpha = I(0.5*noise2))) + geom_line() + ylim(0, max(chains))
g <- g + geom_hline(yintercept = l_proposal)
g
# g <- g + geom_point(data=data.frame(x = hitting_times, y = sapply(1:ncol(chains), function(i) chains[hitting_times[i], i])))

# library(tsibble)
# library(fabletools)
# library(feasts)
# stldf <- chaindf %>% filter(ichain < 500) %>% select(x,y,ichain) %>% group_by(ichain) %>% as_tsibble(index = x, key = ichain) %>%
#     model(STL(y)) %>% components() %>% as_tibble()
# ggplot(stldf, aes(x = x, y = trend, group = ichain)) + geom_line(alpha = 0.1)

##
# xgrid <- seq(from = .1, to = 2, length.out = 250)
# meanhit <- c()
# for (x_ in xgrid){
#     meanhit <- c(meanhit,mean(hitting(l_proposal, x_, nrep = 1000)))
# }
# plot(xgrid, meanhit)
#
# hitting


