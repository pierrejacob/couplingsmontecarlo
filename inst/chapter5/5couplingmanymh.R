## geom ridge lines to show evolution of pair of proposal distributions
## until meeting

rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
set_theme_chapte5 <- function(){
    library(ggplot2)
    library(gridExtra)
    theme_set(theme_void())
    colors <- c(rgb(0.3,0.45,0.8), rgb(0.1, 0.35, 0.35))
    return(list(colors = colors))
}
graphsettings <- set_theme_chapte5()
library(reshape2)
library(dplyr)

rmax_coupling <- function(rp, dp, rq, dq){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)) {
        return(list(xy = c(x, x), identical = TRUE))
    }
    else {
        reject <- TRUE
        y <- NA
        while (reject) {
            y <- rq(1)
            reject <- (dq(y) + log(runif(1)) < dp(y))
        }
        return(list(xy = c(x, y), identical = FALSE))
    }
}
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
    return(rmax_coupling(function(n) rnorm(n, mu1, sigma1),
                         function(x) dnorm(x, mu1, sigma1, log = TRUE),
                         function(n) rnorm(n, mu2, sigma2),
                         function(x) dnorm(x, mu2, sigma2, log = TRUE)))
}

rnorm_max_coupling_givenx <- function(x, mu1, mu2, sigma1, sigma2){
    # rp <- function(n) rnorm(n, mu1, sigma1),
    dp <- function(x) dnorm(x, mu1, sigma1, log = TRUE)
    rq <- function(n) rnorm(n, mu2, sigma2)
    dq <- function(x) dnorm(x, mu2, sigma2, log = TRUE)
    if (dp(x) + log(runif(1)) < dq(x)) {
        return(list(xy = c(x, x), identical = TRUE))
    }
    else {
        reject <- TRUE
        y <- NA
        while (reject) {
            y <- rq(1)
            reject <- (dq(y) + log(runif(1)) < dp(y))
        }
        return(list(xy = c(x, y), identical = FALSE))
    }
}


# target log-pdf
logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
# initial distribution
rinit <- function(){
    chain_state <- rnorm(1, 0, 3)
    current_pdf <- logtarget(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
}
# MH kernel
sd_proposal <- 1
MH_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal <- rnorm(1, chain_state, sd_proposal)
    proposal_pdf <- logtarget(proposal)
    if (log(runif(1)) < (proposal_pdf - current_pdf)){
        return(list(chain_state = proposal, current_pdf = proposal_pdf))
    } else {
        return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
}


## the next function takes a list of states...  WIP
coupledMH_kernel <- function(states){
    # chain_state1 <- state1$chain_state;  current_pdf1 <- state1$current_pdf
    # chain_state2 <- state2$chain_state;  current_pdf2 <- state2$current_pdf
    # proposal from a maximal coupling between successive chains
    nstates <- length(states)
    proposals <- rep(NA, nstates)
    proposal12 <- rnorm_max_coupling(states[[1]]$chain_state, states[[2]]$chain_state, sd_proposal, sd_proposal)
    proposals[1:2] <- proposal12$xy
    sameproposals <- rep(FALSE, nstates-1)
    sameproposals[1] <- proposal12$identical
    if (nstates > 2){
        for (istate in 2:(nstates-1)){
            proposalconsecutive <- rnorm_max_coupling_givenx(proposals[istate], states[[istate]]$chain_state,
                                                             states[[istate+1]]$chain_state, sd_proposal, sd_proposal)
            proposals[istate+1] <- proposalconsecutive$xy[2]
            sameproposals[istate] <- proposalconsecutive$identical
        }
    }
    proposal_pdf <- logtarget(proposals)
    logu <- log(runif(1))
    accepts <- rep(FALSE, nstates)
    for (istate in 1:nstates){
        if (logu < (proposal_pdf[istate] - states[[istate]]$current_pdf)){
            states[[istate]]$current_pdf <- proposal_pdf[istate]
            states[[istate]]$chain_state <- proposals[istate]
            accepts[istate] <- TRUE
        }
    }
    identical_ <- all(sameproposals) && all(accepts)
    return(list(states = states,
                identical = identical_))
}
nstates <- 500
states <- list()
for (istate in 1:nstates) states[[istate]] <- rinit()
nmcmc <- 60
states_history <- matrix(NA, nrow = 1+nmcmc, ncol = nstates)
states_history[1,] <- sapply(states, function(x) x$chain_state)
meetingtime <- Inf
for (imcmc in 1:nmcmc){
    coupled_res <- coupledMH_kernel(states)
    states <- coupled_res$states
    states_history[imcmc+1,] <- sapply(states, function(x) x$chain_state)
    if (is.infinite(meetingtime) && coupled_res$identical) meetingtime <- imcmc
}

states.df <- reshape2::melt(states_history) %>% select(iteration = Var1, ichain = Var2, value = value)
g <- ggplot(data = states.df, aes(x = iteration, y = value, group = ichain)) + geom_line(alpha = 0.2,
                                                                                    colour = graphsettings$colors[1])
g <- g + geom_vline(xintercept = meetingtime+1, colour = graphsettings$colors[2], linetype = 1)
g

ggsave(filename = "../5coupledmh.pdf", plot = g, width = 10, height = 5)

# xgrid <- seq(from = +1, to = 16, length.out = 1e2)
# df1 <- data.frame()
# df2 <- data.frame()
# traces <- data.frame()
# for (time in 1:(coupledchains$meetingtime+1-lag)){
#     heightnorm1 <- sapply(xgrid, function(x) dnorm(x, coupledchains$samples1[time+lag,1], sd_proposal))
#     heightnorm2 <- sapply(xgrid, function(x) dnorm(x, coupledchains$samples2[time,1], sd_proposal))
#     df1 <- rbind(df1, data.frame(x = xgrid, y = time, h = heightnorm1))
#     df2 <- rbind(df2, data.frame(x = xgrid, y = time, h = heightnorm2))
#     traces <- rbind(traces, data.frame(time = time, chain1 = coupledchains$samples1[time+lag,1], chain2 = coupledchains$samples2[time,1]))
# }
# g <- ggplot(df1, aes(x = x, y = y, height = h, group = y)) + geom_ridgeline(alpha = 0.5, fill = graphsettings$colors[1], colour = NA) +
#     geom_ridgeline(data=df2, alpha = 0.5, fill = graphsettings$colors[2], colour = NA) +
#     geom_point(data=traces, aes(group=NULL, x = chain1, y = time, height = NULL), colour = graphsettings$colors[1]) +
#     geom_point(data=traces, aes(group=NULL, x = chain2, y = time, height = NULL), colour = graphsettings$colors[2]) +
#     coord_flip()
# g
#
# # ggsave(plot = g, filename = "../5overlapproposal.pdf", width = 10, height = 5)


