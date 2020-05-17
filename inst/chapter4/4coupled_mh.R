## geom ridge lines to show evolution of pair of proposal distributions
## until meeting

rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter4()
library(ggridges)
library(reshape2)
library(dplyr)
# library(doParallel)
# library(doRNG)
# registerDoParallel(cores = detectCores()-2)

# target log-pdf
logtarget <- function(x) dnorm(x, mean = 0, sd = 10, log = TRUE)
# initial distribution
rinit <- function(){
    chain_state <- rnorm(1, 10, 1)
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

coupledMH_kernel <- function(state1, state2){
    chain_state1 <- state1$chain_state;  current_pdf1 <- state1$current_pdf
    chain_state2 <- state2$chain_state;  current_pdf2 <- state2$current_pdf
    # proposal from a maximal coupling
    proposal <- unbiasedmcmc::rnorm_reflectionmax(chain_state1, chain_state2, sd_proposal)
    proposal_pdf1 <- logtarget(proposal$xy[1])
    # only compute target pdf on 2nd proposal if it is not identical to 1st proposal
    proposal_pdf2 <- proposal_pdf1
    if (!proposal$identical){
        proposal_pdf2 <- logtarget(proposal$xy[2])
    }
    logu <- log(runif(1))
    accept1 <- FALSE; accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
        if (logu < (proposal_pdf1 - current_pdf1)){
            accept1 <- TRUE
            chain_state1 <- proposal$xy[1]; current_pdf1 <- proposal_pdf1
        }
    }
    if (is.finite(proposal_pdf2)){
        if(logu < (proposal_pdf2 - current_pdf2)){
            accept2 <- TRUE
            chain_state2 <- proposal$xy[2]; current_pdf2 <- proposal_pdf2
        }
    }
    identical_ <- proposal$identical && accept1 && accept2
    return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
                state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
                identical = identical_))
}


sample_coupled_chains <- function (single_kernel, coupled_kernel, rinit, m = 1, lag = 1,
                                   max_iterations = Inf, preallocate = 10)
{
    starttime <- Sys.time()
    state1 <- rinit()
    state2 <- rinit()
    dimstate <- length(state1$chain_state)
    nrowsamples1 <- m + preallocate + lag
    samples1 <- matrix(nrow = nrowsamples1, ncol = dimstate)
    samples2 <- matrix(nrow = nrowsamples1 - lag, ncol = dimstate)
    samples1[1, ] <- state1$chain_state
    samples2[1, ] <- state2$chain_state
    time <- 0
    for (t in 1:lag) {
        time <- time + 1
        state1 <- single_kernel(state1)
        samples1[time + 1, ] <- state1$chain_state
    }
    meetingtime <- Inf
    while ((time < max(meetingtime, m)) && (time < max_iterations)) {
        time <- time + 1
        if (is.finite(meetingtime)) {
            state1 <- single_kernel(state1)
            state2 <- state1
        }
        else {
            res_coupled_kernel <- coupled_kernel(state1, state2)
            state1 <- res_coupled_kernel$state1
            state2 <- res_coupled_kernel$state2
            if (res_coupled_kernel$identical) {
                meetingtime <- time
            }
        }
        if ((time + 1) > nrowsamples1) {
            new_rows <- nrowsamples1
            nrowsamples1 <- nrowsamples1 + new_rows
            samples1 <- rbind(samples1, matrix(NA, nrow = new_rows,
                                               ncol = dimstate))
            samples2 <- rbind(samples2, matrix(NA, nrow = new_rows,
                                               ncol = dimstate))
        }
        samples1[time + 1, ] <- state1$chain_state
        samples2[time - lag + 1, ] <- state2$chain_state
    }
    samples1 <- samples1[1:(time + 1), , drop = F]
    samples2 <- samples2[1:(time - lag + 1), , drop = F]
    cost <- lag + 2 * (meetingtime - lag) + max(0, time - meetingtime)
    currenttime <- Sys.time()
    elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) -
                                                         lubridate::ymd_hms(starttime)), "seconds")
    return(list(samples1 = samples1, samples2 = samples2, meetingtime = meetingtime,
                iteration = time, elapsedtime = elapsedtime, cost = cost))
}

lag <- 50
coupledchains <- sample_coupled_chains(MH_kernel, coupledMH_kernel, rinit, lag = lag, m = 200)
coupledchains$meetingtime
g <- qplot(x = 1:nrow(coupledchains$samples1), y = coupledchains$samples1, geom = "blank") + geom_line(colour = graphsettings$colors[1]) +
    geom_line(aes(x = 1:nrow(coupledchains$samples2), y = coupledchains$samples2), colour = graphsettings$colors[2]) +
    geom_vline(xintercept = c(coupledchains$meetingtime, coupledchains$meetingtime+lag), colour = rev(graphsettings$colors))
g

# ggsave(plot = g, filename = "../4laggedcoupled.pdf", width = 10, height = 5)

xgrid <- seq(from = +1, to = 16, length.out = 1e2)
df1 <- data.frame()
df2 <- data.frame()
traces <- data.frame()
for (time in 1:(coupledchains$meetingtime+1-lag)){
    heightnorm1 <- sapply(xgrid, function(x) dnorm(x, coupledchains$samples1[time+lag,1], sd_proposal))
    heightnorm2 <- sapply(xgrid, function(x) dnorm(x, coupledchains$samples2[time,1], sd_proposal))
    df1 <- rbind(df1, data.frame(x = xgrid, y = time, h = heightnorm1))
    df2 <- rbind(df2, data.frame(x = xgrid, y = time, h = heightnorm2))
    traces <- rbind(traces, data.frame(time = time, chain1 = coupledchains$samples1[time+lag,1], chain2 = coupledchains$samples2[time,1]))
}
g <- ggplot(df1, aes(x = x, y = y, height = h, group = y)) + geom_ridgeline(alpha = 0.5, fill = graphsettings$colors[1], colour = NA) +
    geom_ridgeline(data=df2, alpha = 0.5, fill = graphsettings$colors[2], colour = NA) +
    geom_point(data=traces, aes(group=NULL, x = chain1, y = time, height = NULL), colour = graphsettings$colors[1]) +
    geom_point(data=traces, aes(group=NULL, x = chain2, y = time, height = NULL), colour = graphsettings$colors[2]) +
    coord_flip()
g

# ggsave(plot = g, filename = "../4overlapproposal.pdf", width = 10, height = 5)


