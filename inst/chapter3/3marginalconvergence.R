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


# logtarget <- function(x) log(0.5 * dnorm(x, mean = 0, sd = 1) + 0.5 * dnorm(x, mean = -4, sd = 1))
# logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
logtarget <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = c(1.2,2.1), log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
curve(sapply(x, function(xp) exp(logtarget(xp))), from = -5, to = 8)
sd_proposal <- .5
## density p
rinit <- function(){
  if (runif(1) < 0.3){
    rnorm(1,10,0.5)
  } else {
    rnorm(1,-5,3)
  }
}

# xsamples <- sapply(1:1e3, function(index) rinit())
# hist(xsamples, nclass = 100)

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

nchains <- 5000
niterations <- 150
# chains <- matrix(0, nrow = niterations, ncol = nchains)
chains <- foreach(ichain = 1:nchains, .combine = cbind) %dorng% {
  chain <- rep(0, niterations)
  chain_state <- rinit()
  for (i in 1:niterations){
    chain_state <- MH_kernel(chain_state)
    chain[i] <- chain_state
  }
  chain
}

# matplot(chains, type = "l", col = "black")
# hist(chains[niterations,], nclass = 100)

chains <- matrix(chains, nrow = niterations)
chains.df <- melt(chains)
names(chains.df) <- c("iteration", "chain", "value")
chains.df %>% tail
g <- ggplot(chains.df, aes(x = value, y = factor(iteration), height = ..density..)) +
  geom_density_ridges(scale = 10, fill = graphsettings$colors[2], stat = "density", adjust = .5)
g <- g + scale_y_discrete(breaks = c(0,50,100,150, 200))
g <- g + xlab("x") + ylab("iteration") + xlim(-14, 14)
g <- g + coord_flip()
g
ggsave(filename = "../mcmc.convergence.ridge.pdf", plot = g, width = 10, height = 7)

