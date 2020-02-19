rm(list = ls())
source("commongrounds.R")
## snippets of code 
## compute average of weights on log-scale
N <- 1e4
x <- rnorm(N)
logf <- function(x) dcauchy(x, log = TRUE)
logg <- function(x) dnorm(x, log = TRUE)
logw <- logf(x) - logg(x) 
mlw <- max(logw) 
mlw + log(mean(exp(logw - mlw)))

## with some probability do this, otherwise do that
if (log(runif(1)) < (logf(x[1]) - logg(x[1]))){
  Z <- -Inf 
} else {
  Z <- x
}


## sample a bivariate Markov chain
T <- 50
x <- matrix(0, nrow = T, ncol = 2)
xcurrent <- c(1,-1)
for (t in 1:T){ 
  increment <- rnorm(1)
  xcurrent <- 0.9 * xcurrent + increment
  x[t,] <- xcurrent
} 

## inverse cdf transform
nsamples <- 1e4
u <- runif(nsamples, min = 0, max = 1)
x <- -log(u)
##

gscatter <- qplot(x = x, y = u, geom = "blank") + geom_point(alpha = 0.05) 
gmargx <- qplot(x = x, geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[2]) + stat_function(fun = dexp)
gmargy <- qplot(x = u, geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[1]) + stat_function(fun = dunif) + coord_flip()
empty <- ggplot()
g <- gridExtra::grid.arrange(gmargx, empty, gscatter, gmargy, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
# ggsave(filename = "icdf.pdf", plot = g)

## rejection sampling
nsamples <- 1e3
x <- mvtnorm::rmvnorm(nsamples, mean = c(0,0), sigma = matrix(c(1,0.8,0.8,1), nrow = 2, ncol = 2))
g <- qplot(x[,1], x[,2], geom = "blank") + geom_point() + xlim(-4, 4) + ylim(-4,4)
# constraint: x + y <= 1
rtrunc <- function(){
  accepted <- FALSE
  while (!accepted){
    x <- mvtnorm::rmvnorm(1, mean = c(0,0), sigma = matrix(c(1,0.8,0.8,1), nrow = 2, ncol = 2))
    if (x[1,1] + x[1,2] <= 1){
      return(x)
    }
  }
}
for (isample in 1:nsamples) x[isample,] <- rtrunc()
g <- qplot(x[,1], x[,2], geom = "blank") + geom_point() + xlim(-4, 4) + ylim(-4,4)
g <- g + geom_abline(slope = -1, intercept = 1)
g
## or...
## sample from a beta(3,2) by sampling from uniforms

minimizer <- optimize(f = function(x) -dbeta(x, 3, 2, log = T), interval = c(0,1))$minimum
logM <- dbeta(minimizer, 3, 2, log = T)
#
y <- runif(nsamples)
us <- runif(nsamples)
accepts <- log(us) < (dbeta(y, 3, 2, log = TRUE) - 0 - logM)
x <- y[accepts]

pt.df <- data.frame(u = us, y = y, accept = accepts)
g <- ggplot(pt.df, aes(x = y, y = us  * exp(logM), colour = accept))  + geom_point() + 
  scale_color_manual(values = colors) + stat_function(fun = function(x) dbeta(x, 3, 2), colour = "black") +
  theme(legend.position = "none")
g
# ggsave(filename = "rejection.pdf", plot = g)

## importance sampling
rq <- function(n) rgamma(n, 2, 1)
logdq <- function(x) dgamma(x, 2, 1, log = TRUE)
logdp <- function(x) dgamma(x, 4, 2, log = TRUE)
nsamples <- 1e4
x <- rq(nsamples)
logw <- logdp(x) - logdq(x)
mlw <- max(logw) 
mlw + log(mean(exp(logw - mlw)))
w <- exp(logw - mlw)
nw <- w / sum(w) # normalized weights

ghist1 <- qplot(x = x, geom = "blank") + geom_histogram(aes(y=..density..), alpha = 1, fill = colors[1]) +
  stat_function(fun = function(x) exp(logdq(x)), colour = colors[1]) + ylim(0,.6) + xlim(0,12)
  # geom_histogram(aes(y=..density.., weight = nw), alpha = 0.75, fill = colors[2]) +
  # stat_function(fun = function(x) exp(logdp(x)), colour = colors[2]) 

glogw <- qplot(x = x, xend = x, y = 0, yend = logw, geom = "blank") + geom_segment(alpha = 1) + xlim(0,12) + ylim(-7,2)

ghist2 <- qplot(x = x, geom = "blank") + geom_histogram(aes(y=..density.., weight = nw), alpha = 1, fill = colors[2]) +
  stat_function(fun = function(x) exp(logdp(x)), colour = colors[2]) + ylim(0,.6) + xlim(0,12)

g <- gridExtra::grid.arrange(ghist1, glogw, ghist2, ncol = 1)
# ggsave(filename = "importancesampling.pdf", plot = g)



## Metropolis--Hastings
## target: Normal with mean c(0,0) and with variance (1, 0.8, 0.8, 1)
## transition kernel
dimension <- 2
logdtarget <- function(x) mvtnorm::dmvnorm(x, mean = c(0,0), sigma = matrix(c(1,0.8,0.8,1), nrow = 2), log = TRUE)
nmcmc <- 10000
## initialization
state <- list(x = c(5, -5))
state$logdtarget <- logdtarget(state$x)
## store entire history
chain <- matrix(nrow = dimension, ncol = nmcmc)
proposal_history <- matrix(nrow = dimension, ncol = nmcmc)
## number of accepted proposals
naccepts <- 0
for (imcmc in 1:nmcmc){
  ## propose new point
  proposal <- state$x + rnorm(dimension, mean = 0, sd = 1) 
  proposal_history[,imcmc] <- proposal
  ## compute target density at proposed point
  proposal_logdtarget <- logdtarget(proposal)
  ## decide to accept or not
  accept <- log(runif(1)) < proposal_logdtarget - state$logdtarget
  if (accept){
    ## if accept
    naccepts <- naccepts + 1
    state <- list(x = proposal, logdtarget = proposal_logdtarget)
  }
  ## else, no modification to state
  ## store current state
  chain[,imcmc] <- state$x
}
cat("accept rate", naccepts / nmcmc, "%\n")
gtrace1 <- qplot(x = 1:100, y = chain[1,1:100], geom = "blank") + geom_line()
gtrace1 <- gtrace1 + geom_point(aes(y = proposal_history[1,1:100]), colour = colors[1]) + geom_point(aes(y = chain[1,1:100]), colour = colors[2])

ghist1 <- qplot(x = chain[1,100:nmcmc], geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[2]) + stat_function(fun = dnorm)

gtrace2 <- qplot(x = 1:100, y = chain[2,1:100], geom = "blank") + geom_line()
gtrace2 <- gtrace2 + geom_point(aes(y = proposal_history[2,1:100]), colour = colors[1]) + geom_point(aes(y = chain[2,1:100]), colour = colors[2])
ghist2 <- qplot(x = chain[2,100:nmcmc], geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[2]) + stat_function(fun = dnorm)
g <- gridExtra::grid.arrange(gtrace1, ghist1, gtrace2, ghist2, ncol=2, nrow=2, widths=c(4, 1), heights=c(2, 2))
# ggsave(filename = "metropolis.pdf", plot = g)

## Gibbs sampling on Ising model

cppFunction('
  IntegerMatrix ising_gibbs_sweep_(IntegerMatrix state, NumericVector proba_beta){
  RNGScope scope;
  int size = state.rows();
  int s;
  int itop, ibottom, jright, jleft;
  GetRNGstate();
  for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
      s = 0;
      itop = (i+1) % size;
      ibottom = ((i + size - 1) % size);
      jright = (j+1) % size;
      jleft = (j + size - 1) % size;
      s += state(itop, j) + state(ibottom, j) + state(i, jright) + state(i, jleft);
      state(i,j) = 2*((runif(1))(0) < proba_beta((s+4)/2)) - 1;
    }
  }
  PutRNGstate();
  return state;
}
')

beta <- 0.42
ss <- c(-4, -2, 0, 2, 4)
n <- 30
df <- expand.grid(1:n, 1:n)
probas <- exp(ss * beta) / (exp(ss * beta) + exp(-ss * beta))
state <- matrix(2*(runif(n*n) < 0.5)-1, nrow = n)
image(state)
df$z <- as.numeric(state)
g1 <- ggplot(df, aes(x = Var1, y = Var2, fill = factor(z))) + geom_tile() + scale_fill_manual(values = c("white", "pink")) +
  theme(legend.position = "none")

for (imcmc in 1:2) state <- ising_gibbs_sweep_(state, probas)
df$z <- as.numeric(state)
g2 <-ggplot(df, aes(x = Var1, y = Var2, fill = factor(z))) + geom_tile() + scale_fill_manual(values = c("white", "pink")) +
  theme(legend.position = "none")

for (imcmc in 1:1e2) state <- ising_gibbs_sweep_(state, probas)
df$z <- as.numeric(state)
g3 <- ggplot(df, aes(x = Var1, y = Var2, fill = factor(z))) + geom_tile() + scale_fill_manual(values = c("white", "pink")) +
  theme(legend.position = "none")

g <- gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
# ggsave(filename = "gibbs.pdf", plot = g, width = 15, height = 5)

