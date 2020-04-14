rm(list=ls())

library(couplingsmontecarlo)
library(numDeriv)
graphsettings <- set_theme_chapter1()

## Metropolis--Hastings
## target: Normal with mean c(0,0) and with variance (1, 0.8, 0.8, 1)
## transition kernel
dimension <- 2
logdtarget <- function(x) mvtnorm::dmvnorm(x,
                                          mean=c(0,0),
                                          sigma=matrix(c(1, 0.8, 0.8, 1),
                                                       nrow=2),
                                          log=TRUE)
nmcmc <- 10000
## initialization
state <- list(x=c(5, -5))
state$logdtarget <- logdtarget(state$x)
## store entire history
chain <- matrix(nrow=dimension, ncol=nmcmc)
proposal_history <- matrix(nrow=dimension, ncol=nmcmc)
## number of accepted proposals
naccepts <- 0
for (imcmc in 1:nmcmc){
    ## propose new point
    proposal <- state$x + rnorm(dimension, mean=0, sd=1)
    proposal_history[,imcmc] <- proposal
    ## compute target density at proposed point
    proposal_logdtarget <- logdtarget(proposal)
    ## decide to accept or not
    accept <- log(runif(1)) < proposal_logdtarget - state$logdtarget
    if (accept){
        ## if accept
        naccepts <- naccepts + 1
        state <- list(x=proposal, logdtarget=proposal_logdtarget)
    }
    ## else, no modification to state
    ## store current state
    chain[, imcmc] <- state$x
}
cat("accept rate", naccepts / nmcmc, "%\n")
gtrace1 <- qplot(x=1:100, y=chain[1, 1:100], geom="blank") +
            geom_line() +
            geom_point(aes(y=proposal_history[1, 1:100]), colour=graphsettings$colors[1]) +
            geom_point(aes(y=chain[1, 1:100]), colour=graphsettings$colors[2])
ghist1 <- qplot(x=chain[1, 100:nmcmc], geom="blank") +
            geom_histogram(aes(y=..density..), fill=graphsettings$colors[2]) +
            stat_function(fun=dnorm)
gtrace2 <- qplot(x=1:100, y=chain[2, 1:100], geom="blank") +
            geom_line() +
            geom_point(aes(y=proposal_history[2, 1:100]), colour=graphsettings$colors[1]) +
            geom_point(aes(y=chain[2, 1:100]), colour=graphsettings$colors[2])
ghist2 <- qplot(x=chain[2, 100:nmcmc], geom="blank") +
            geom_histogram(aes(y=..density..), fill=graphsettings$colors[2]) +
            stat_function(fun=dnorm)
gridExtra::grid.arrange(gtrace1, ghist1, gtrace2, ghist2,
                            ncol=2, nrow=2,
                            widths=c(4, 1), heights=c(2, 2))
# ggsave(filename="metropolis.pdf", plot=g)

## Gibbs sampling on Ising model

src <-
    "IntegerMatrix ising_gibbs_sweep_(IntegerMatrix state, NumericVector proba_beta){
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
                s += state(itop, j)
                    + state(ibottom, j)
                    + state(i, jright)
                    + state(i, jleft);
                state(i,j) = 2*((runif(1))(0) < proba_beta((s+4)/2)) - 1;
            }
        }
        PutRNGstate();
        return state;
    }"
Rcpp::cppFunction(src)

beta <- 0.42
ss <- c(-4, -2, 0, 2, 4)
n <- 30
df <- expand.grid(1:n, 1:n)
probas <- exp(ss * beta) / (exp(ss * beta) + exp(-ss * beta))
state <- matrix(2 * (runif(n * n) < 0.5) - 1, nrow=n)
image(state)
df$z <- as.numeric(state)
g1 <- ggplot(df, aes(x=Var1, y=Var2, fill=factor(z))) +
        geom_tile() +
        scale_fill_manual(values=c("white", "pink")) +
        theme(legend.position="none")

for (imcmc in 1:2) state <- ising_gibbs_sweep_(state, probas)
df$z <- as.numeric(state)
g2 <- ggplot(df, aes(x=Var1, y=Var2, fill=factor(z))) +
        geom_tile() +
        scale_fill_manual(values=c("white", "pink")) +
        theme(legend.position="none")

for (imcmc in 1:1e2) state <- ising_gibbs_sweep_(state, probas)
df$z <- as.numeric(state)
g3 <- ggplot(df, aes(x=Var1, y=Var2, fill=factor(z))) +
        geom_tile() +
        scale_fill_manual(values=c("white", "pink")) +
        theme(legend.position="none")

gridExtra::grid.arrange(g1, g2, g3, nrow=1)
# ggsave(filename="gibbs.pdf", plot=g, width=15, height=5)
