rm(list = ls())
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter1()
set.seed(1)

## we follow quite closely the nice blog post written here
## https://freakonometrics.hypotheses.org/60575
## by Arthur Charpentier (https://freakonometrics.github.io)

## the blog post provides data
# loc_fichier = "http://freakonometrics.free.fr/titanic.RData"
# download.file(loc_fichier, "titanic.RData")
# load("titanic.RData")

## alternately get data from the following package, source being Kaggle.com
# install.packages("titanic")
library(titanic)
dim(titanic_train)
dim(titanic_test)
head(titanic::titanic_train)
head(titanic::titanic_test)
# ??titanic
## let's define our data using titanic_train, remove rows with missing data, removing the column "name"
summary(titanic_train)
library(dplyr)
df <- na.omit(titanic_train)
X <- df %>% select(Pclass, Sex, Age, SibSp)
Y <- df$Survived

mean(Y)
mean(Y[X$Pclass==1])
mean(Y[X$Pclass==2])
mean(Y[X$Pclass==3])



table(X$Sex)
X$Sex <- sapply(X$Sex, function(x) ifelse(x == "male", 1, 0))
X$Pclass2 <- sapply(X$Pclass, function(x) ifelse(x == "2", 1, 0))
X$Pclass3 <- sapply(X$Pclass, function(x) ifelse(x == "3", 1, 0))
X <- X %>% select(-Pclass)
X$Age <- (X$Age - mean(X$Age))/(sd(X$Age))
# X$Age2 <- X$Age^2
# X$Age3 <- X$Age^3
X <- X %>% select(Pclass2, Pclass3, Sex, Age, SibSp)
head(X)
glm_regression <- glm(Y ~ X$Pclass2 + X$Pclass3 + X$Sex + X$Age + X$SibSp, family = "binomial", data = NULL)
summary(glm_regression)

## manually add column for the intercept
X <- cbind(1, X)

##
# mean(X$Sex[X$Pclass3==0])
# mean(X$Sex[X$Pclass3==1])

library(Rcpp)
## function to compute the log-likelihood and its gradient, for each 'beta' stored as a column of 'betas'
## Y is a vector of zeros and ones, and X a matrix of covariates with nrow = length(Y)
## returns list with 'gradients' and 'lls'
cppFunction('
List loglikelihood(const Eigen::MatrixXd & betas, const Eigen::ArrayXd & Y, const Eigen::MatrixXd & X){
  int p = X.cols();
  int n = X.rows();
  int nbetas = betas.cols();
  double eval = 0;
  Eigen::MatrixXd xbeta = X * betas;
  Eigen::ArrayXXd expxbeta = xbeta.array().exp();
  Eigen::ArrayXXd gradients(p,nbetas);
  gradients.setZero(p,nbetas);
  Eigen::ArrayXd lls(nbetas);
  lls.setZero();
  for (int i = 0; i < n; i ++){
    lls += Y(i) * xbeta.row(i).array() - (1 + expxbeta.row(i)).log();
    for (int j = 0; j < nbetas; j ++){
      gradients.col(j) += X.row(i).array() * (Y(i) - expxbeta(i,j) / (1. + expxbeta(i,j)));
    }
  }
  return List::create(Named("gradients") = wrap(gradients), Named("evals") = wrap(lls));
}', depends = "RcppEigen")

## Normal prior with same variance sigma2 on each coefficient
## 'betas' is a matrix with p rows and 'nbetas' columns
## returns gradients (matrix p x nbetas) and evaluation of prior log-density
logprior <- function(betas, sigma2){
  lps <- - 0.5 * colSums(betas * betas) / sigma2 - 0.5 * dim(betas)[1] * log(sigma2)
  gradients <- - betas / sigma2
  return(list(gradients = gradients, evals = lps))
}

# we will set the prior variance to
sigma2prior <- 5^2
## MLE
beta0 <- as.numeric(glm_regression$coefficients) + rnorm(length(glm_regression$coefficients), sd = 0.001)
# loglikelihood(matrix(beta0, ncol = 1), Y, as.matrix(X))
loglikehood_optim_results <- optim(par = beta0, function(beta) -loglikelihood(matrix(beta, ncol = 1), Y, as.matrix(X))$evals)
## compare parameters obtained by glm or 'manual' optimization of the likelihood
loglikehood_optim_results$par
glm_regression$coefficients
## pretty close match
## compute Hessian at MLE
hessian_at_mle <- numDeriv::hessian(func = function(beta) -loglikelihood(matrix(beta, ncol = 1), Y, as.matrix(X))$evals, loglikehood_optim_results$par)
## thus Laplace approximation of the posterior has mean and variance
post_mean_laplace <- loglikehood_optim_results$par
post_cov_laplace <- solve(hessian_at_mle)

## importance sampling from Laplace approximation
nsamples_is <- 1e4
samples_is <- t(mvtnorm::rmvnorm(nsamples_is, post_mean_laplace, post_cov_laplace))
samples_is_ll <- loglikelihood(samples_is, Y, as.matrix(X))
samples_is_lp <- logprior(samples_is, sigma2prior)$evals
logweights <- samples_is_ll$evals + samples_is_lp - mvtnorm::dmvnorm(t(samples_is), post_mean_laplace, post_cov_laplace, log=TRUE)
# hist(logweights)
mlw <- max(logweights)
avew <- mlw + log(mean(exp(logweights - mlw))) # weight average
w <- exp(logweights - mlw) # weights
nw <- w / sum(w) # normalized weights
ess <- 1/sum(nw^2)
cat("Effective sample size:", floor(ess), "out of", nsamples_is, "draws")
## seems to be a very good importance sampling proposal
post_mean_is <- colSums(t(samples_is) * nw)

## we can also plot a marginal distribution
## qplot(x = samples_is[3,], weight = nw, geom = "blank") + geom_histogram() + theme_minimal()

## initialize chains
init_chains <- function(nchains){
  # samples_ <- t(mvtnorm::rmvnorm(nchains, post_mean_laplace, post_cov_laplace))
  samples_ <- t(mvtnorm::rmvnorm(nchains, rep(0, ncol(X)), diag(1, ncol(X), ncol(X))))
  samples_ll <- loglikelihood(samples_, Y, as.matrix(X))
  samples_lp <- logprior(samples_, sigma2prior)
  return(list(states=samples_, lls = samples_ll, lps = samples_lp))
}
## the 'chains" contain "states" (actual values of the regression coefficients)
## as well as evaluated log-likelihood, gradients thereof, log prior and gradients thereof

## function to perform one step of MH for each of the chains
rwmh_kernel <- function(chains, proposal_variance){
  # number of chains
  nchains <- ncol(chains$states)
  # propose new values
  proposals <- chains$states + t(mvtnorm::rmvnorm(nchains, rep(0, nrow(chains$states)), proposal_variance))
  # compute log-likelihood and log-prior
  proposal_lls <- loglikelihood(proposals, Y, as.matrix(X))$evals
  proposal_lps <- logprior(proposals, sigma2prior)$evals
  # compute acceptance ratios
  logacceptratios <- (proposal_lls + proposal_lps) - (chains$lls$evals + chains$lps$evals)
  # draw uniforms and compare them to acceptance ratios
  accepts <- (log(runif(nchains)) < logacceptratios)
  for (ichain in 1:nchains){
    if (accepts[ichain]){
      # replace current states by proposed ones
      chains$states[,ichain] <- proposals[,ichain]
      chains$lls$evals[ichain] <- proposal_lls[ichain]
      chains$lps$evals[ichain] <- proposal_lps[ichain]
    }
  }
  return(chains)
}

# run eight chains
nchains <- 8
# initialize
chains_current <- init_chains(nchains)
# we won't need the gradients for MH
chains_current$lls$gradients <- NULL
chains_current$lps$gradients <- NULL
# try certain proposal covariance matrix
proposal_variance <- post_cov_laplace
# total number of iterations
niterations_mh <- 5000
# store all states of the chains
chains_history_mh <- array(dim = c(1+niterations_mh, nrow(chains_current$states), nchains))
chains_history_mh[1,,] <- chains_current$states
for (iteration in 1:niterations_mh){
  chains_current <- rwmh_kernel(chains_current, proposal_variance)
  chains_history_mh[iteration+1,,] <- chains_current$states
}
## trace plot of all the chains, first component, first 1000 iterations
chains_beginning <- chains_history_mh[1:1000,1,]
chains_beginning <- reshape2::melt(chains_beginning)
##
gtraceplot_mh <- ggplot(chains_beginning, aes(x = Var1, group = Var2, y = value)) + geom_line(size = 0.3)
ggsave(filename = "../logistitanic-traceplot-mh.pdf", plot = gtraceplot_mh, width = 10, height = 5)
## choice of burn-in, arbitrary
burnin_mh <- 1000
## seems stationary, visually
## hacky way of retrieving average acceptance rate
acceptrate <- 1 - mean(abs(diff(chains_history_mh[(burnin_mh+1):(niterations_mh+1),1,]))<1e-10)
cat("MH acceptance rate", acceptrate*100, "%\n")
## acceptance rate is neither too close to zero or one, so can be deemed alright

## compare posterior obtained with IS and with random walk MH
imarg <- 3
qplot(x = samples_is[imarg,], weight = nw, geom = "blank") + geom_density() + theme_minimal() +
  geom_density(aes(x = as.numeric(chains_history_mh[(burnin_mh+1):(niterations_mh+1),imarg,]), weight = NULL), linetype = 2) +
  xlab(paste0("coefficient ", imarg))
## very close agreement of the methods

## stack all chains together
allchains_mh <- do.call(rbind, lapply(1:nchains, function(ichain) chains_history_mh[(burnin_mh+1):(niterations_mh+1),,ichain]))
## compare posterior mean estimates
colMeans(allchains_mh)
post_mean_is
## very close agreement

## create a "pair plot"
library(GGally)
my_fn <- function(data, mapping, ...){
  print(mapping)
  p <- ggplot(data = data, mapping = mapping) +
    stat_density_2d(contour = TRUE, colour = "black") +
    scale_x_continuous(breaks = seq(from = -5, to = 5, by = 0.1))
    # stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
    # viridis::scale_fill_viridis()
  p
}

cor_fun <- function(data, mapping, method="pearson", ndp=2, sz=5, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- 3+sz* abs(est)
  lbl <- round(est, ndp)
  ggplot(data=data, mapping=mapping) +
    annotate("text", x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE), label=lbl, size=lb.size,...)+
    theme(panel.grid = element_blank())
}

## takes a while to run
g <- ggpairs(data.frame(allchains_mh), lower=list(continuous=my_fn),
        diag=list(continuous=wrap("barDiag")),
        upper=list(continuous=wrap(cor_fun, sz=5, stars=FALSE)),
        axisLabels = "show",
        columnLabels = c("intercept", "class 2", "class 3", "Sex", "Age", "# Siblings" )) +
        theme(strip.text.x = element_text(size = 12, vjust = 1),
              axis.text.x = element_text(size = 5, angle = -45))
# g
# ggsave(filename = "../logistitanic-pairs.pdf", plot = g, width = 10, height = 10)

## next we can predict the survival probabilities of the two
## main characters in Titanic

newX = data.frame(
  Pclass = as.factor(c(1,3)),
  Sex = as.factor(c("female","male")),
  Age = c(17,20),
  SibSp = c(1,0))
## Rose is female, 17 years old, first class passenger, travels with her mother
## Jack is male, 20, 3rd class, travels with no siblings

newX$Sex <- sapply(newX$Sex, function(x) ifelse(x == "male", 1, 0))
newX$Pclass2 <- sapply(newX$Pclass, function(x) ifelse(x == "2", 1, 0))
newX$Pclass3 <- sapply(newX$Pclass, function(x) ifelse(x == "3", 1, 0))
newX <- newX %>% select(-Pclass)
## note: rescaling of Age must be identical to that of training data
newX$Age <- (newX$Age - mean(df$Age))/(sd(df$Age))
# newX$Age2 <- newX$Age^2
# newX$Age3 <- newX$Age^3
newX <- newX %>% select(Pclass2, Pclass3, Sex, Age, SibSp)
newX <- cbind(1, as.matrix(newX))

## compute 'beta X' for each beta produced by the MH chain
beta_times_newX <- allchains_mh %*% t(newX)
## transform into probability of survival
survival_prob_newX <- colMeans(exp(beta_times_newX)/(1+exp(beta_times_newX)))
cat("survival probabilities for Rose and Jack:", 100 * survival_prob_newX, "% respectively\n")

## function to perform HMC moves on each chain
hmc_kernel <- function(chains, leapfrognsteps, leapfrogepsilon, massmatrix){
  nchains <- ncol(chains$states)
  dimstate <- nrow(chains$states)
  ## generate momenta variables
  initial_momenta <- t(mvtnorm:::rmvnorm(nchains, rep(0, dimstate), massmatrix))
  grad_ <- chains$lls$gradients + chains$lps$gradients
  positions <- chains$states
  trajectories <- array(dim = c(leapfrognsteps, dim(positions)))
  ## leap frog integrator
  momenta <- initial_momenta + leapfrogepsilon * grad_ / 2
  for (step in 1:leapfrognsteps){
    positions <- positions + leapfrogepsilon * massmatrix_inv %*% momenta
    eval_prior <- logprior(positions, sigma2 = sigma2prior)
    eval_ll    <- loglikelihood(positions, Y, as.matrix(X))
    if (step != leapfrognsteps){
      momenta <- momenta + leapfrogepsilon * (eval_prior$gradients + eval_ll$gradients)
    }
    trajectories[step,,] <- positions
  }
  momenta <- momenta + leapfrogepsilon * (eval_prior$gradients + eval_ll$gradients) / 2
  ## Now MH acceptance step
  proposed_pdfs <- eval_prior$evals     +  eval_ll$evals
  current_pdfs  <- chains$lps$evals + chains$lls$evals
  mhratios <- proposed_pdfs - current_pdfs
  mhratios <- mhratios + (-0.5 * colSums(momenta * (massmatrix_inv %*% momenta))) - (-0.5 * colSums(initial_momenta * (massmatrix_inv %*% initial_momenta)))
  if (any(is.na(mhratios))) mhratios[is.na(mhratios)] <- -Inf
  accepts <- log(runif(nchains)) < mhratios
  accept_rate <- mean(accepts)
  for (ichain in 1:nchains){
    if (accepts[ichain]){
      chains$states[,ichain]              <- positions[,ichain]
      chains$lps$evals[ichain]    <- eval_prior$evals[ichain]
      chains$lps$gradients[,ichain] <- eval_prior$gradients[,ichain]
      chains$lls$evals[ichain]              <- eval_ll$evals[ichain]
      chains$lls$gradients[,ichain]    <- eval_ll$gradients[,ichain]
    }
  }
  return(list(chains = chains, accepts = accepts, trajectories = trajectories))
}

nchains <- 6
chains_current <- init_chains(nchains)
leapfrognsteps <- 10
leapfrogepsilon <- 0.1
massmatrix <- solve(post_cov_laplace)
massmatrix_inv <- post_cov_laplace
niterations_hmc <- 2000
chains_history_hmc <- array(dim = c(1+niterations_hmc, nrow(chains_current$states), nchains))
chains_history_hmc[1,,] <- chains_current$states
alltrajectories_firstchain <- matrix(nrow = niterations_hmc * leapfrognsteps, ncol = nrow(chains_current$states))
for (iteration in 1:niterations_hmc){
  hmc_results <- hmc_kernel(chains_current, leapfrognsteps, leapfrogepsilon, massmatrix)
  chains_current <- hmc_results$chains
  chains_history_hmc[iteration+1,,] <- chains_current$states
  alltrajectories_firstchain[((iteration-1)*leapfrognsteps+1):(iteration*leapfrognsteps),] <- hmc_results$trajectories[,,1]
}

## traceplot
chains_beginning <- chains_history_hmc[1:250,1,]
chains_beginning <- reshape2::melt(chains_beginning)
##
gtraceplot_hmc <- ggplot(chains_beginning, aes(x = Var1, group = Var2, y = value)) + geom_line(size = 0.3)
gtraceplot_hmc

## burnin
burnin_hmc <- 100
acceptrate <- 1-mean(abs(diff(chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),1,]))<1e-10)
cat("HMC acceptance rate", acceptrate*100, "%\n")


imarg <- 1
qplot(x = samples_is[imarg,], weight = nw, geom = "blank") + geom_density() + theme_minimal() +
  geom_density(aes(x = as.numeric(chains_history_mh[(burnin_mh+1):(niterations_mh+1),imarg,]), weight = NULL), linetype = 2) +
  xlab(paste0("coefficient ", imarg)) +
  geom_density(aes(x = as.numeric(chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),imarg,]), weight = NULL), linetype = 3)


apply(chains_history_mh[(burnin_mh+1):(niterations_mh+1),,], 2, mean)
apply(chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),,], 2, mean)
post_mean_is

apply(chains_history_mh[(burnin_mh+1):(niterations_mh+1),,], 2, sd)
apply(chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),,], 2, sd)

## visualize full trajectory with leap frog intermediate steps
alltraj_df <- data.frame(alltrajectories_firstchain)
alltraj_df$irow <- 1:nrow(alltraj_df)
alltraj_df <- alltraj_df %>% mutate(iteration = floor(irow/leapfrognsteps))
head(alltraj_df, 12)

allchains_hmc <- do.call(rbind, lapply(1:nchains, function(ichain) chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),,ichain]))
dim(allchains_hmc)

gtraj <- ggplot(data = data.frame(allchains_hmc), aes(x = X2, y = X3)) +
  stat_density_2d(contour = TRUE, colour =  graphsettings$colors[1])
gtraj <- gtraj + geom_path(data = alltraj_df %>% filter(irow >= 100, irow <= 400), colour = graphsettings$colors[2])
# gtraj <- gtraj + geom_path(data = alltraj_df %>% filter(irow > 100, irow < 400), arrow = arrow(length = unit(.1, "inch")), aes(group = iteration))
gtraj <- gtraj + geom_point(data = alltraj_df %>% filter(irow >= 100, irow <= 400, irow %% leapfrognsteps == 0), colour = graphsettings$colors[2])
gtraj
ggsave(filename = "../logistitanic-hmctraj.pdf", plot = gtraj, width = 7, height = 7)

## comparison of autocorrelograms

maxlag <- 100
lags <- 0:maxlag
acfvalues_mh <- acf(chains_history_mh[(burnin_mh+1):(niterations_mh+1),1,1], plot = FALSE, lag.max = maxlag)$acf
acfvalues_hmc <- acf(chains_history_hmc[(burnin_hmc+1):(niterations_hmc+1),1,1], plot = FALSE, lag.max = maxlag)$acf
g <- qplot(x = lags, y = acfvalues_mh, geom = "blank")
g <- g + geom_segment(aes(xend = lags, yend = 0), col = graphsettings$colors[1])
g <- g + geom_segment(aes(x = lags+0.5, xend = lags+0.5, y = acfvalues_hmc, yend = 0), col = graphsettings$colors[2])
g <- g + ylab("ACF") + xlab("Lag")
g

ggsave(filename = "../logistitanic-acfcomparison.pdf", plot = g, width = 7, height = 7)


# ## impact of burn-in?
# ## no burnin
# chain1_mh <- chains_history_mh[1:(niterations_mh+1-burnin_mh),,1]
# # with burnin
# chain1_mh_burnt <- chains_history_mh[(burnin_mh+1):(niterations_mh+1),,1]
# matplot(apply(chain1_mh, 2, cumsum) / 1:nrow(chain1_mh), type = 'l', col = 'black')
# matplot(apply(chain1_mh_burnt, 2, cumsum) / 1:nrow(chain1_mh_burnt), type = 'l', col = 'red', add = TRUE)



