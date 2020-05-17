
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
logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
# initial distribution
rinit <- function(){
    chain_state <- rnorm(1, 0, 10)
    current_pdf <- logtarget(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
}
sd_proposal <- .5
# MH kernel
coupledMH_kernel <- function(state1, state2, coupledproposal){
    chain_state1 <- state1$chain_state;  current_pdf1 <- state1$current_pdf
    chain_state2 <- state2$chain_state;  current_pdf2 <- state2$current_pdf
    # proposal from a maximal coupling
    proposal <- coupledproposal(chain_state1, chain_state2, sd_proposal)
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

run_coupled_chains <- function(niterations, coupledproposal){
    chain1 <- rinit()
    chain2 <- rinit()
    chain1_history <- rep(NA, niterations+1)
    chain1_history[1] <- chain1$chain_state
    chain2_history <- rep(NA, niterations+1)
    chain2_history[1] <- chain2$chain_state
    for (iteration in 1:niterations){
        c_results <- coupledMH_kernel(chain1, chain2, coupledproposal)
        chain1 <- c_results$state1
        chain2 <- c_results$state2
        chain1_history[iteration+1] <- chain1$chain_state
        chain2_history[iteration+1] <- chain2$chain_state
    }
    return(list(chain1_history = chain1_history,
                chain2_history = chain2_history))
}

niterations <- 200
indep_coupledproposal <- function(state1, state2, sd_proposal){
    list(xy = c(rnorm(1, mean = state1, sd = sd_proposal), rnorm(1, mean = state2, sd = sd_proposal)), identical = FALSE)
}
crn_coupledproposal <- function(state1, state2, sd_proposal){
    increment <- rnorm(1, mean = 0, sd = sd_proposal)
    list(xy = c(state1 + increment, state2 + increment), identical = FALSE)
}
max_coupledproposal <- function(state1, state2, sd_proposal){
    unbiasedmcmc::rnorm_max_coupling(state1, state2, sd_proposal, sd_proposal)
}
reflmax_coupledproposal <- function(state1, state2, sd_proposal){
    unbiasedmcmc::rnorm_reflectionmax(state1, state2, sd_proposal)
}

## repeat a number of times, and plot squared distances between chains
nrep <- 1000
indep_cc <- foreach(irep = 1:nrep) %dorng% {
    run_coupled_chains(niterations, indep_coupledproposal)
}
crn_cc <- foreach(irep = 1:nrep) %dorng% {
    run_coupled_chains(niterations, crn_coupledproposal)
}
max_cc <- foreach(irep = 1:nrep) %dorng% {
    run_coupled_chains(niterations, max_coupledproposal)
}
reflmax_cc <- foreach(irep = 1:nrep) %dorng% {
    run_coupled_chains(niterations, reflmax_coupledproposal)
}

indep_cc.df <- reshape2::melt(do.call(what = rbind, args = lapply(indep_cc, function(l) abs(l$chain1_history - l$chain2_history)))) %>%
    rename(ichain = Var1, iteration = Var2, distance = value) %>% mutate(coupling = "indep")
crn_cc.df <- reshape2::melt(do.call(what = rbind, args = lapply(crn_cc, function(l) abs(l$chain1_history - l$chain2_history)))) %>%
    rename(ichain = Var1, iteration = Var2, distance = value) %>% mutate(coupling = "crn")
max_cc.df <- reshape2::melt(do.call(what = rbind, args = lapply(max_cc, function(l) abs(l$chain1_history - l$chain2_history)))) %>%
    rename(ichain = Var1, iteration = Var2, distance = value) %>% mutate(coupling = "reflection")
reflmax_cc.df <- reshape2::melt(do.call(what = rbind, args = lapply(reflmax_cc, function(l) abs(l$chain1_history - l$chain2_history)))) %>%
    rename(ichain = Var1, iteration = Var2, distance = value) %>% mutate(coupling = "reflection max")

g <- ggplot(indep_cc.df %>% group_by(iteration) %>% summarise(meandist = mean(distance)), aes(x = iteration, y = meandist)) +
    geom_line(colour = graphsettings$colors[1]) +
    geom_line(data = reflmax_cc.df %>% group_by(iteration) %>% summarise(meandist = mean(distance)), colour = graphsettings$colors[1], linetype = 2) +
    geom_line(data = crn_cc.df %>% group_by(iteration) %>% summarise(meandist = mean(distance)), colour = graphsettings$colors[2]) +
    geom_line(data = max_cc.df %>% group_by(iteration) %>% summarise(meandist = mean(distance)), colour = graphsettings$colors[2], linetype = 2) +
    scale_y_log10()
# ggsave(filename = "../4distancebetweenchains.diffcouplings.pdf", plot = g, width = 10, height = 5)
