## remove objects from current memory
## and load packages
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
## package for the titanic data set
library(titanic)
## let's define our data using titanic_train, remove rows with missing data, removing the column "name"
summary(titanic_train)
df <- na.omit(titanic_train)
X <- df %>% select(Pclass, Sex, Age, SibSp)
Y <- df$Survived

table(X$Sex)
X$Sex <- sapply(X$Sex, function(x) ifelse(x == "male", 1, 0))
X$Pclass2 <- sapply(X$Pclass, function(x) ifelse(x == "2", 1, 0))
X$Pclass3 <- sapply(X$Pclass, function(x) ifelse(x == "3", 1, 0))
X <- X %>% select(-Pclass)
X$Age <- (X$Age - mean(X$Age))/(sd(X$Age))
# X$Age2 <- X$Age^2
# X$Age3 <- X$Age^3
X <- X %>% select(Pclass2, Pclass3, Sex, Age, SibSp)
glm_regression <- glm(Y ~ X$Pclass2 + X$Pclass3 + X$Sex + X$Age + X$SibSp, family = "binomial", data = NULL)
summary(glm_regression)

## manually add column for the intercept
X <- cbind(1, X)

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

##### next, MCMC part
## initial distribution of the chains
init_chains <- function(nchains){
    # samples_ <- t(mvtnorm::rmvnorm(nchains, post_mean_laplace, post_cov_laplace))
    samples_ <- t(mvtnorm::rmvnorm(nchains, rep(0, ncol(X)), diag(1, ncol(X), ncol(X))))
    samples_ll <- loglikelihood(samples_, Y, as.matrix(X))
    samples_lp <- logprior(samples_, sigma2prior)
    return(list(states=samples_, lls = samples_ll, lps = samples_lp))
}
## the 'chains" contain "states" (actual values of the regression coefficients)
## as well as evaluated log-likelihood, gradients thereof, log prior and gradients thereof
## HMC kernel performs one step of Hamiltonian Monte Carlo
hmc_kernel <- function(chain, tuning){
    attach(tuning, warn.conflicts = F)
    dimstate <- nrow(chain$states)
    ## generate momenta variables
    initial_momenta <- t(mvtnorm:::rmvnorm(1, rep(0, dimstate), massmatrix))
    grad_ <- chain$lls$gradients + chain$lps$gradients
    positions <- chain$states
    ## leap frog integrator
    momenta <- initial_momenta + leapfrogepsilon * grad_ / 2
    for (step in 1:leapfrognsteps){
        positions <- positions + leapfrogepsilon * massmatrix_inv %*% momenta
        eval_prior <- logprior(positions, sigma2 = sigma2prior)
        eval_ll    <- loglikelihood(positions, Y, as.matrix(X))
        if (step != leapfrognsteps){
            momenta <- momenta + leapfrogepsilon * (eval_prior$gradients + eval_ll$gradients)
        }
    }
    momenta <- momenta + leapfrogepsilon * (eval_prior$gradients + eval_ll$gradients) / 2
    ## Now MH acceptance step
    proposed_pdfs <- eval_prior$evals     +  eval_ll$evals
    current_pdfs  <- chain$lps$evals + chain$lls$evals
    mhratios <- proposed_pdfs - current_pdfs
    mhratios <- mhratios + (-0.5 * colSums(momenta * (massmatrix_inv %*% momenta))) - (-0.5 * colSums(initial_momenta * (massmatrix_inv %*% initial_momenta)))
    if (any(is.na(mhratios))) mhratios[is.na(mhratios)] <- -Inf
    accept <- log(runif(1)) < mhratios
    if (accept){
        chain$states              <- positions
        chain$lps$evals    <- eval_prior$evals
        chain$lps$gradients <- eval_prior$gradients
        chain$lls$evals             <- eval_ll$evals
        chain$lls$gradients    <- eval_ll$gradients
    }
    return(chain)
}
##
## function to perform HMC moves on two chains, using same noise
coupled_hmc_kernel <- function(chain1, chain2, tuning){
    attach(tuning, warn.conflicts = F)
    nchains <- 2
    dimstate <- nrow(chain1$states)
    ## generate common momenta variables
    initial_momenta <- t(mvtnorm:::rmvnorm(1, rep(0, dimstate), massmatrix))
    ##
    grad1_ <- chain1$lls$gradients + chain1$lps$gradients
    grad2_ <- chain2$lls$gradients + chain2$lps$gradients
    position1 <- chain1$states
    position2 <- chain2$states
    ## leap frog integrator
    momenta1 <- initial_momenta + leapfrogepsilon * grad1_ / 2
    momenta2 <- initial_momenta + leapfrogepsilon * grad2_ / 2
    for (step in 1:leapfrognsteps){
        position1 <- position1 + leapfrogepsilon * massmatrix_inv %*% momenta1
        position2 <- position2 + leapfrogepsilon * massmatrix_inv %*% momenta2
        eval_prior <- logprior(cbind(position1, position2), sigma2 = sigma2prior)
        eval_ll    <- loglikelihood(cbind(position1, position2), Y, as.matrix(X))
        if (step != leapfrognsteps){
            momenta1 <- momenta1 + leapfrogepsilon * (eval_prior$gradients[,1] + eval_ll$gradients[,1])
            momenta2 <- momenta2 + leapfrogepsilon * (eval_prior$gradients[,2] + eval_ll$gradients[,2])
        }
    }
    momenta1 <- momenta1 + leapfrogepsilon * (eval_prior$gradients[,1] + eval_ll$gradients[,1]) / 2
    momenta2 <- momenta2 + leapfrogepsilon * (eval_prior$gradients[,2] + eval_ll$gradients[,2]) / 2
    ## Now MH acceptance step for each chain
    proposed_pdfs <- eval_prior$evals     +  eval_ll$evals
    current_pdfs  <- c(chain1$lps$evals, chain2$lps$evals) + c(chain1$lls$evals, chain2$lls$evals)
    u_ar <- runif(1)
    ## first chain
    mhratio <- proposed_pdfs[1] - current_pdfs[1]
    mhratio <- mhratio + (-0.5 * colSums(momenta1 * (massmatrix_inv %*% momenta1))) - (-0.5 * colSums(initial_momenta * (massmatrix_inv %*% initial_momenta)))
    if (is.na(mhratio)) mhratio <- -Inf
    if (log(u_ar) < mhratio){
        chain1$states <- position1
        chain1$lls$evals <- eval_ll$evals[1]
        chain1$lls$gradients <- eval_ll$gradients[,1,drop=F]
        chain1$lps$evals <- eval_prior$evals[1]
        chain1$lps$gradients <- eval_prior$gradients[,1,drop=F]
    }
    ## second chain
    mhratio <- proposed_pdfs[2] - current_pdfs[2]
    mhratio <- mhratio + (-0.5 * colSums(momenta2 * (massmatrix_inv %*% momenta2))) - (-0.5 * colSums(initial_momenta * (massmatrix_inv %*% initial_momenta)))
    if (is.na(mhratio)) mhratio <- -Inf
    if (log(u_ar) < mhratio){
        chain2$states <- position2
        chain2$lls$evals <- eval_ll$evals[2]
        chain2$lls$gradients <- eval_ll$gradients[,2,drop=F]
        chain2$lps$evals <- eval_prior$evals[2]
        chain2$lps$gradients <- eval_prior$gradients[,2,drop=F]
    }
    return(list(chain1 = chain1, chain2 = chain2, identical = FALSE))
}
## some tuning parameters for HMC
## number of leap frog steps
leapfrognsteps <- 10
## step size of the leap frog integrator
leapfrogepsilon <- 0.1
## mass matrix = covariance of the momentum variable
## taken as approximation of the posterior precision matrix
massmatrix <- solve(post_cov_laplace)
## inverse of mass matrix
massmatrix_inv <- post_cov_laplace
tuning <- list(leapfrognsteps = leapfrognsteps, leapfrogepsilon = leapfrogepsilon,
               massmatrix = massmatrix, massmatrix_inv = massmatrix_inv)

## run pair of HMC chains propagated with same noise
run_coupled_hmc <- function(niterations, tuning){
    chain1 <- init_chains(1)
    chain1 <- hmc_kernel(chain1, tuning)
    chain2 <- init_chains(1)
    chain1_history <- matrix(NA, nrow = niterations+1, ncol = nrow(chain1$states))
    chain1_history[1,] <- chain1$states
    chain2_history <- matrix(NA, nrow = niterations, ncol = nrow(chain2$states))
    for (iteration in 1:niterations){
        chmc_results <- coupled_hmc_kernel(chain1, chain2, tuning)
        chain1 <- chmc_results$chain1
        chain2 <- chmc_results$chain2
        chain1_history[iteration+1,] <- chain1$states
        chain2_history[iteration,] <- chain2$states
    }
    return(list(chain1_history = chain1_history,
                chain2_history = chain2_history))
}
niterations <- 30
chains_history_hmc <- run_coupled_hmc(niterations, tuning)
chain1_history.df <- reshape2::melt(chains_history_hmc$chain1_history) %>% select(iteration=Var1, component=Var2, value) %>%
    mutate(ichain = 1, iteration = iteration - 1) %>% filter(iteration > 0)
chain2_history.df <- reshape2::melt(chains_history_hmc$chain2_history) %>% select(iteration=Var1, component=Var2, value) %>%
    mutate(ichain = 2)
## we can see chains contracting to one another very rapidly
ggplot(data = rbind(chain1_history.df, chain2_history.df),
       aes(x = iteration, y = value, colour = factor(ichain), group = interaction(component, ichain))) + geom_line() +
    theme_minimal() + theme(legend.position = "none")  +
    scale_color_manual(values = graphsettings$colors)

## repeat a number of times, and plot squared distances between chains
nrep <- 100
indep_coupled_chains <- foreach(irep = 1:nrep) %dorng% {
    run_coupled_hmc(niterations, tuning)
}
## compute squared euclidean distance between chains
squared_distances <- lapply(indep_coupled_chains, function(coupledchains){
    squareddist_ <- apply(coupledchains$chain1_history[2:nrow(coupledchains$chain1_history),] - coupledchains$chain2_history, 1, function(v) sum(v^2))
})
distances.df <- reshape2::melt(matrix(unlist(squared_distances), ncol = nrep)) %>% rename(iteration = Var1, ichain = Var2, distance = value)
## plot squared distances against iterations
ggplot(distances.df,
       aes(x = iteration, y = distance, group = ichain)) + geom_line(alpha = 0.2) + theme_minimal() +
    ylab("squared distance") +
    scale_y_log10()
## so after 30 steps or so pairs of chains tend to be very close to one another

distsummary <- distances.df %>% group_by(iteration) %>% summarise(logd = log(mean(distance)))
ggplot(distsummary, aes(x = iteration, y = logd)) + geom_line() + theme_bw()
logd_slope <- lm(logd ~ iteration, data = distsummary %>% as.data.frame)$coef['iteration']
logd_slope
## logd_t = C + logd_slope * t
## so ... d_t = Ctilde * exp(logd_slope * t)
geom_rate <- exp(logd_slope)
print(geom_rate)

## Rhee Glynn estimator
rtruncation <- function(pgeom) 1 + rgeom(1, pgeom)
truncationtailproba <- function(x, pgeom) 1-pgeom(x-2, pgeom)
# pgeom <- 0.123
# xx <- sapply(1:1e5, function(l) rtruncation(pgeom))
# mean(xx>=0)
# truncationtailproba(0, pgeom)
# mean(xx>=1)
# truncationtailproba(1, pgeom)
# mean(xx>=2)
# truncationtailproba(2, pgeom)
# mean(xx>=3)
# truncationtailproba(3, pgeom)
# mean(xx>=4)
# truncationtailproba(4, pgeom)


rheeglynn <- function(pgeom){
    truncation <- rtruncation(pgeom)
    cc <- run_coupled_hmc(niterations = truncation, tuning)
    delta <- matrix(NA, nrow = truncation+1, ncol = ncol(cc$chain1_history))
    delta[1,] <- cc$chain1_history[1,]
    for (t in 1:truncation){
        delta[t+1,] <- cc$chain1_history[t+1,] - cc$chain2_history[t,]
        delta[t+1,] <- delta[t+1,] / truncationtailproba(t, pgeom)
    }
    return(list(truncation = truncation, estimator = colSums(delta)))
}


nsamples <- 1e3

pgeom <- 0.95
rg_estimators <- foreach(irep = 1:nsamples) %dorng% {
    rheeglynn(pgeom)
}

rg_estimators_matrix <- lapply(rg_estimators, function(l) data.frame(matrix(l$estimator, nrow = 1))) %>% bind_rows()
hist(rg_estimators_matrix[,1], nclass = 100)
mean(rg_estimators_matrix[,1])
var(rg_estimators_matrix[,1])

# pgeom_seq <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
pgeom_seq <- c(0.6, 0.8)
rgmeancost <- rep(0, length(pgeom_seq))
rgvar <- matrix(0, nrow = 6, ncol = length(pgeom_seq))
rgestimator <- matrix(0, nrow = 6, ncol = length(pgeom_seq))

for (ip in seq_along(pgeom_seq)){
    print(ip)
    pgeom <- pgeom_seq[ip]
    rg_estimators <- foreach(irep = 1:nsamples) %dorng% {
        rheeglynn(pgeom)
    }
    rg_estimators_matrix <- lapply(rg_estimators, function(l) data.frame(matrix(l$estimator, nrow = 1))) %>% bind_rows()
    rgestimator[,ip] <- colMeans(rg_estimators_matrix)
    rgvar[,ip] <- apply(rg_estimators_matrix, 2, var) / nsamples
    rgmeancost[ip] <- mean(1 + 2*sapply(rg_estimators, function(l) l$truncation))
}

# plot(pgeom_seq, rgmeancost, type = 'l')
plot(pgeom_seq, rgvar[1,], type = 'l')
# plot(pgeom_seq, rgmeancost * rgvar[1,], type = 'l')
# plot(pgeom_seq, rgestimator[1,], type = 'l')




