## this implements a coupling
## of the sampler in Section 6.2 of Rapidly Mixing Markov Chains: A Comparison of Techniques  by Venkatesan Guruswami
## sampling uniformly subsets of size k from set of n elements
rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
set_theme_chapter4 <- function(){
    library(ggplot2)
    library(gridExtra)
    theme_set(theme_void())
    colors <- c(rgb(0.8,0.5,0.2), rgb(0.2, 0.6, 0.9))
    return(list(colors = colors))
}
graphsettings <- set_theme_chapter4()
library(ggridges)
library(reshape2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
## number of elements in total
n <- 10
fullset <- 1:n
## number of elements in subset
k <- 3
## initial state
rinit <- function() 1:k
state_current <- rinit()
## Markov kernel
single_kernel <- function(state){
    u <- runif(1)
    if (u < 0.5){
        return(state)
    } else {
        i <- sample(x = state, size = 1)
        j <- sample(x = setdiff(fullset, state), size = 1)
        return(c(setdiff(state, i), j))
    }
}
## MCMC run for 'nmcmc' iterations
nmcmc <- 1e4
state_history <- matrix(NA, nrow = nmcmc, ncol = k)
state_history[1,] <- state_current
for (imcmc in 2:nmcmc){
    state_current <- single_kernel(state_current)
    state_history[imcmc,] <- state_current
}
##
## trace plot ...
df <- data.frame(iteration = 1:50)
df$label <- apply(state_history[df$iteration,], 1, function(v) paste0("{", paste(sort(v)-1, collapse = ","), "}"))
ggplot(df, aes(y = iteration, x = 0, label = label)) + geom_text() + scale_y_reverse()

## Now, how do we know whether the chains converged?
burnin <- 1e3
table(as.numeric(state_history[burnin:nmcmc,])) / ((nmcmc-burnin+1)*k)
## the frequencies of occurrence of each element look pretty close to uniform ('1/n')
## for small n and k we can actually list all the possible subsets
## number of subsets with k elements out of n
choose(n, k)
all_subsets <- combn(1:n, k, simplify = F)
## all_subsets is a list of 'choose(n, k)' vectors of sorted indices

sorted_state_history <- t(apply(state_history, 1, sort))
freq_subsets <- rep(0, choose(n,k))
freq_subsets <- foreach (isubset = 1:choose(n,k), .combine = c) %dopar% {
    subset <- all_subsets[[isubset]]
    ## count how many times this subset appeared post-burnin
    freq_subset_ <- sum(apply(sorted_state_history[burnin:nmcmc,], 1, function(v) all(v == subset)))
    freq_subset_ / (nmcmc-burnin+1)
}
plot(freq_subsets, ylim = c(0, 2/choose(n,k)))
abline(h = 1/choose(n,k))
## at this point the sampler seems to be doing OK but we still don't know how much were actually
## necessary for burn-in

## function to sample two elements from maximal coupling of uniform
## distributions on subsets of {1,...,n}
max_unif_sets <- function(subset1, subset2){
    s1 <- length(subset1)
    s2 <- length(subset2)
    # original vectors of probabilities
    p1 <- rep(0, n)
    p2 <- rep(0, n)
    p1[subset1] <- 1/s1
    p2[subset2] <- 1/s2
    # decompose into common part
    common_part <- pmin(p1, p2)
    c <- sum(common_part)
    # and residual parts
    r1 <- (p1-common_part)/(1-c)
    r2 <- (p2-common_part)/(1-c)
    # sample pair of indices
    if (runif(1) < c){
        index <- sample(x = 1:n, size = 1, prob = common_part/c)
        return(c(index, index))
    } else {
        index1 <- sample(x = 1:n, size = 1, prob = r1)
        index2 <- sample(x = 1:n, size = 1, prob = r2)
        return(c(index1, index2))
    }
}

## coupled kernel
coupled_kernel <- function(state1, state2){
    u <- runif(1)
    if (u < 0.5){
        return(list(state1 = state1, state2 = state2, identical = FALSE))
    } else {
        ## sample i and iprime from maximal coupling of uniform distributions
        i12 <- max_unif_sets(state1, state2)
        i1 <- i12[1]; i2 <- i12[2]
        ## sample i and iprime from maximal coupling of uniform distributions
        j12 <- max_unif_sets(setdiff(fullset, state1), setdiff(fullset, state2))
        j1 <- j12[1]; j2 <- j12[2]
        return(list(state1 = c(setdiff(state1, i1), j1),
                    state2 = c(setdiff(state2, i2), j2)))
    }
}

## construction of two chains with a lag L
## using single_kernel and coupled_kernel
sample_meetingtime <- function(single_kernel, coupled_kernel, rinit, lag = 1, max_iterations = Inf){
    starttime <- Sys.time()
    # initialize two chains
    state1 <- rinit(); state2 <- rinit()
    # move first chain for 'lag' iterations
    time <- 0
    for (t in 1:lag){
        time <- time + 1
        state1 <- single_kernel(state1)
    }
    # move two chains until meeting (or until max_iterations)
    meetingtime <- Inf
    # two chains could already be identical by chance
    if (all(sort(state1) == sort(state2))) meetingtime <- lag
    while (is.infinite(meetingtime) && (time < max_iterations)){
        time <- time + 1
        # use coupled kernel
        coupledstates <- coupled_kernel(state1, state2)
        state1 <- coupledstates$state1
        state2 <- coupledstates$state2
        # check if meeting has occurred
        if (all(sort(state1) == sort(state2))) meetingtime <- time
    }
    currentime <- Sys.time()
    elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
    return(list(meetingtime = meetingtime, elapsedtime = elapsedtime))
}

## generate a number of meeting times, for a certain lag
nrep <- 1e3
lag <- 30
meetings <- foreach(irep = 1:nrep) %dorng% sample_meetingtime(single_kernel, coupled_kernel, rinit, lag)
meeting_times <- sapply(meetings, function(x) x$meetingtime)
## plot meeting times
hist(meeting_times - lag)
## compute TV upper bounds
tv_upper_bound_estimates <- function(meeting_times, L, t){
    return(mean(pmax(0,ceiling((meeting_times-L-t)/L))))
}
niter <- 50
upperbounds <- sapply(1:niter, function(t) tv_upper_bound_estimates(unlist(meeting_times), lag, t))
g_tvbounds <- qplot(x = 1:niter, y = upperbounds, geom = "line")
g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
g_tvbounds <- g_tvbounds + labs(title = "uniform sampling of subsets of size k") + ylim(0,1)
g_tvbounds + theme_minimal()

## this suggest that 20 or 30 steps might be enough
## we can check that by producing independent chains for this many steps
## and check whether we have a good approximation of the posterior
nrep <- 1e4
nmcmc <- 20
parallel_chains <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    state_current <- rinit()
    for (imcmc in 1:nmcmc){
        state_current <- single_kernel(state_current)
    }
    state_current
}
## frequencies of each element
table(as.numeric(parallel_chains)) / (nrep*k)
## compute frequencies of each subset
sorted_parallel_chains <- t(apply(parallel_chains, 1, sort))
freq_subsets <- rep(0, choose(n,k))
freq_subsets <- foreach (isubset = 1:choose(n,k), .combine = c) %dopar% {
    subset <- all_subsets[[isubset]]
    ## count how many times this subset appeared post-burnin
    freq_subset_ <- sum(apply(sorted_parallel_chains, 1, function(v) all(v == subset)))
    freq_subset_ / nrep
}
plot(freq_subsets, ylim = c(0, 2/choose(n,k)))
abline(h = 1/choose(n,k))



