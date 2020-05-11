## this implements a coupling
## of the sampler in Section 6.2 of Rapidly Mixing Markov Chains: A Comparison of Techniques  by Venkatesan Guruswami
## sampling uniformly subsets of size k from set of n elements
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
## number of elements in total
n <- 10
fullset <- 0:(n-1)
## number of elements in subset
k <- 3
## initial state
rinit <- function() 0:(k-1)
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

state_current <- rinit()
print(state_current)

state_current <- single_kernel(state_current)
print(state_current)

## MCMC run for 'nmcmc' iterations
nmcmc <- 1e3
## and 'nchains' independent chains
nchains <- 8
chains <- foreach(ichain = 1:nchains) %dorng% {
    state_history <- matrix(NA, nrow = nmcmc, ncol = k)
    state_current <- rinit()
    state_history[1,] <- state_current
    for (imcmc in 2:nmcmc){
        state_current <- single_kernel(state_current)
        state_history[imcmc,] <- state_current
    }
    state_history
}

chaindf <- lapply(seq_along(chains), function(ichain){
    data.frame(ichain = ichain,
               iteration = 1:nrow(chains[[ichain]]),
               label = apply(chains[[ichain]], 1, function(v) paste0("{", paste(sort(v), collapse = ","), "}")))
})
chaindf <- do.call(rbind, chaindf)

tail(chaindf)
gtraceplot <- ggplot(chaindf %>% filter(iteration < 30),
                     aes(y = iteration, x = ichain, label = label, group = ichain,
                         colour = factor(ichain %% 2)))
gtraceplot <- gtraceplot + geom_text() + scale_y_reverse()
gtraceplot <- gtraceplot + scale_color_manual(values = graphsettings$colors)
gtraceplot <- gtraceplot + theme(legend.position = "none")
gtraceplot
# ggsave(filename = "../traceplot.subsets.pdf", plot = gtraceplot)

## Now, how do we know whether the chains have converged?
## for small n and k we can actually list all the possible subsets
## number of subsets with k elements out of n
choose(n, k)
all_subsets <- combn(fullset, k, simplify = F)
## all_subsets is a list of 'choose(n, k)' vectors of sorted indices
all_subsets <- lapply(all_subsets, function(v) paste0("{", paste(sort(v), collapse = ","), "}"))
burnin <- 1e2


freq_subsets <- rep(0, choose(n,k))
freq_subsets <- foreach (isubset = 1:choose(n,k), .combine = c) %dopar% {
    subset <- all_subsets[[isubset]]
    ## count how many times this subset appeared post-burnin
    count_subset <- sum(chaindf %>% filter(iteration > burnin) %>% mutate(count_subset = (label == subset)) %>% pull(count_subset))
    count_subset / ((nmcmc-burnin+1) * nchains)
}
plot(freq_subsets, ylim = c(0, 2/choose(n,k)))
abline(h = 1/choose(n,k))
## at this point the sampler seems to be doing OK but we still don't know how much were actually
## necessary for burn-in

## we can check that by producing independent chains for this many steps
## and check whether we have a good approximation of the posterior
nchains <- 1e4
nmcmc <- 10
parallel_chains <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
    state_current <- rinit()
    for (imcmc in 1:nmcmc){
        state_current <- single_kernel(state_current)
    }
    c(ichain, paste0("{", paste(sort(state_current), collapse = ","), "}"))
}
parallel_chains <- data.frame(ichain = parallel_chains[,1], label = parallel_chains[,2], row.names = NULL)
head(parallel_chains)
## compute frequencies of each subset
freq_subsets <- rep(0, choose(n,k))
freq_subsets <- foreach (isubset = 1:choose(n,k), .combine = c) %dopar% {
    subset <- all_subsets[[isubset]]
    ## count how many times this subset appeared post-burnin
    count_subset <- sum(parallel_chains %>% mutate(count_subset = (label == subset)) %>% pull(count_subset))
    count_subset / nchains
}
plot(freq_subsets, ylim = c(0, 2/choose(n,k)))
abline(h = 1/choose(n,k))



