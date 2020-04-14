rm(list = ls())
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter1()
set.seed(1)
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

Rcpp::cppFunction(
"int ising_sum_states(const IntegerMatrix & state){
  int size = state.rows();
  int s = 0;
  int i1;
  int j1;
  for (int i = 0; i < size; i++){
    if (i == (size - 1)){
      i1 = 0;
    } else {
      i1 = i+1;
    }
    for (int j = 0; j < size; j++){
      if (j == (size - 1)){
        j1=0;
      } else {
        j1 = j+1;
      }
      s += state(i,j) * (state(i,j1) + state(i1,j));
    }
  }
  return s;
}")

beta <- 0.4
ss <- c(-4, -2, 0, 2, 4)
n <- 30
df <- expand.grid(1:n, 1:n)
probas <- exp(ss * beta) / (exp(ss * beta) + exp(-ss * beta))
state <- matrix(2 * (runif(n * n) < 0.5) - 1, nrow=n)
df$z <- as.numeric(state)
g1 <- ggplot(df, aes(x=Var1, y=Var2, fill=factor(z))) +
    geom_tile() +
    scale_fill_manual(values=c(graphsettings$colors[1], graphsettings$colors[2])) +
    theme(legend.position="none")
g1
niterations <- 3e4
state_history <- array(dim = c(niterations, dim(state)))
for (imcmc in 1:niterations){
    state <- ising_gibbs_sweep_(state, probas)
    state_history[imcmc,,] <- state
}


df$z <- as.numeric(state)
g2 <- ggplot(df, aes(x=Var1, y=Var2, fill=factor(z))) +
    geom_tile() +
    scale_fill_manual(values=c(graphsettings$colors[1], graphsettings$colors[2])) +
    theme(legend.position="none")

g <- gridExtra::grid.arrange(g1, g2, nrow=1)
ggsave(filename="../isingstates.pdf", plot=g, width=10, height=5)

plot(sapply(1:1e3, function(index) ising_sum_states(state_history[index,,])), type = 'l')
burnin <- 1e3

proba_11 <- cumsum(state_history[(burnin+1):(niterations),1,1])
proba_11 <- proba_11 / (1:length(proba_11))
matplot(proba_11, type = 'l')

