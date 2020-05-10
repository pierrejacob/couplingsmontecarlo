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



