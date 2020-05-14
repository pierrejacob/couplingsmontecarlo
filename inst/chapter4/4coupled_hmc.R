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
    chain2 <- init_chains(1)
    chain1_history <- matrix(NA, niterations, nrow(chain1$states))
    chain2_history <- matrix(NA, niterations, nrow(chain2$states))
    for (iteration in 1:niterations){
        chmc_results <- coupled_hmc_kernel(chain1, chain2, tuning)
        chain1 <- chmc_results$chain1
        chain2 <- chmc_results$chain2
        chain1_history[iteration,] <- chain1$states
        chain2_history[iteration,] <- chain2$states
    }
    return(list(chain1_history = chain1_history,
                chain2_history = chain2_history))
}
niterations <- 30
chains_history_hmc <- run_coupled_hmc(niterations, tuning)
chain1_history.df <- reshape2::melt(chains_history_hmc$chain1_history) %>% select(iteration=Var1, component=Var2, value) %>% mutate(ichain = 1)
chain2_history.df <- reshape2::melt(chains_history_hmc$chain2_history) %>% select(iteration=Var1, component=Var2, value) %>% mutate(ichain = 2)
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
    squareddist_ <- apply(coupledchains$chain1_history - coupledchains$chain2_history, 1, function(v) sum(v^2))
})
distances.df <- reshape2::melt(matrix(unlist(squared_distances), ncol = nrep)) %>% rename(iteration = Var1, ichain = Var2, distance = value)
## plot squared distances against iterations
ggplot(distances.df,
       aes(x = iteration, y = distance, group = ichain)) + geom_line(alpha = 0.2) + theme_minimal() +
  ylab("squared distance") +
  scale_y_log10()
## so after 30 steps or so pairs of chains tend to be very close to one another

## Next we construct a mixture of kernels so that chains can exactly meet
## the mixture has one component equal to the HMC kernel implemented above
## the other component is a random walk MH kernel with Normal proposals
## because we can easily couple such RWMH kernel to generate exact meetings,
## via maximal couplings... specifically, 'reflection-maximal' couplings here

## function to perform MH step for one chain
rwmh_kernel <- function(chain, tuning){
    # propose new values
    proposal <- mvtnorm::rmvnorm(1, mean = chain$states[,1], tuning$Sigma)
    # compute log-likelihood and log-prior
    proposal_lls <- loglikelihood(t(proposal), Y, as.matrix(X))
    proposal_lps <- logprior(t(proposal), sigma2prior)
    # compute acceptance ratios
    logacceptratio <- (proposal_lls$evals + proposal_lps$evals) - (chain$lls$evals[1] + chain$lps$evals[1])
    # draw uniforms and compare them to acceptance ratios
    accept <- (log(runif(1)) < logacceptratio)
    if (accept){
        # replace current states by proposed ones
        chain$states[,1] <- proposal
        chain$lls$evals[1] <- proposal_lls$evals
        chain$lls$gradients[,1] <- proposal_lls$gradients
        chain$lps$evals[1] <- proposal_lps$evals
        chain$lps$gradients[,1] <- proposal_lps$gradients
    }
    return(chain)
}
## function to perform coupled MH step for two chains
## with 'reflection-maximal' coupling of the proposals
coupled_rwmh_kernel <- function(chain1, chain2, tuning){
    # propose new values
    proposals <- couplingsmontecarlo:::rmvnorm_reflection_max_coupling_(chain1$states[,1], chain2$states[,1],
                                                                        tuning$Sigma_chol, tuning$inv_Sigma_chol)
    proposals_states <- proposals$xy
    # compute log-likelihood and log-prior
    proposal_lls <- loglikelihood(proposals_states, Y, as.matrix(X))
    proposal_lps <- logprior(proposals_states, sigma2prior)
    # compute acceptance ratios
    logacceptratios <- (proposal_lls$evals + proposal_lps$evals) -
        (c(chain1$lls$evals, chain2$lls$evals) + c(chain1$lps$evals, chain2$lps$evals))
    # draw uniforms and compare them to acceptance ratios
    accepts <- (log(runif(2)) < logacceptratios)
    if (accepts[1]){
        chain1$states[,1] <- proposals_states[,1,drop=F]
        chain1$lls$evals[1] <- proposal_lls$evals[1]
        chain1$lls$gradients[,1] <- proposal_lls$gradients[,1,drop=F]
        chain1$lps$evals[1] <- proposal_lps$evals[1]
        chain1$lps$gradients[,1] <- proposal_lps$gradients[,1,drop=F]
    }
    if (accepts[2]){
        chain2$states[,1] <- proposals_states[,2,drop=F]
        chain2$lls$evals[1] <- proposal_lls$evals[2]
        chain2$lls$gradients[,1] <- proposal_lls$gradients[,2,drop=F]
        chain2$lps$evals[1] <- proposal_lps$evals[2]
        chain2$lps$gradients[,1] <- proposal_lps$gradients[,2,drop=F]
    }
    i_ <- (proposals$identical && all(accepts))
    return(list(chain1 = chain1, chain2 = chain2, identical = i_))
}
##
## add tuning parameters relating to covariance of Normal proposals
tuning$Sigma <- diag(rep(1e-6, ncol(X)))
tuning$Sigma_chol <- chol(tuning$Sigma)
tuning$inv_Sigma_chol <- solve(chol(tuning$Sigma))

## mixture of the kernels
tuning$hmc_weight <- 0.8
## single kernel
mixkernel <- function(chain, tuning){
    if (runif(1) < tuning$hmc_weight){
        return(hmc_kernel(chain, tuning))
    } else {
        return(rwmh_kernel(chain, tuning))
    }
}
## coupled kernel
coupled_mixkernel <- function(chain1, chain2, tuning){
  if (runif(1) < tuning$hmc_weight){
    return(coupled_hmc_kernel(chain1, chain2, tuning))
  } else {
    return(coupled_rwmh_kernel(chain1, chain2, tuning))
  }
}

## function to sample coupled chains until they meet
sample_coupled_chains <- function(single_kernel, coupled_kernel, rinit, tuning, m = 1, lag = 1,
                                   max_iterations = Inf, preallocate = 10){
  starttime <- Sys.time()
  state1 <- rinit(1)
  state2 <- rinit(1)
  dimstate <- nrow(state1$states)
  nrowsamples1 <- m + preallocate + lag
  samples1 <- matrix(nrow = nrowsamples1, ncol = dimstate)
  samples2 <- matrix(nrow = nrowsamples1 - lag, ncol = dimstate)
  samples1[1, ] <- state1$states
  samples2[1, ] <- state2$states
  time <- 0
  for (t in 1:lag) {
    time <- time + 1
    state1 <- single_kernel(state1, tuning)
    samples1[time + 1, ] <- state1$states
  }
  meetingtime <- Inf
  while ((time < max(meetingtime, m)) && (time < max_iterations)) {
    time <- time + 1
    if (is.finite(meetingtime)) {
      state1 <- single_kernel(state1, tuning)
      state2 <- state1
    }
    else {
      res_coupled_kernel <- coupled_kernel(state1, state2, tuning)
      state1 <- res_coupled_kernel$chain1
      state2 <- res_coupled_kernel$chain2
      if (res_coupled_kernel$identical) {
        meetingtime <- time
      }
    }
    if ((time + 1) > nrowsamples1) {
      new_rows <- nrowsamples1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows,
                                         ncol = dimstate))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows,
                                         ncol = dimstate))
    }
    samples1[time + 1, ] <- state1$states
    samples2[time - lag + 1, ] <- state2$states
  }
  samples1 <- samples1[1:(time + 1), , drop = F]
  samples2 <- samples2[1:(time - lag + 1), , drop = F]
  cost <- lag + 2 * (meetingtime - lag) + max(0, time - meetingtime)
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) -
                                                     lubridate::ymd_hms(starttime)), "seconds")
  return(list(samples1 = samples1, samples2 = samples2, meetingtime = meetingtime,
              iteration = time, elapsedtime = elapsedtime, cost = cost))
}

## run pairs of chains until they meet
## a number of times, in parallel
coupled_chains_tilmeet <- foreach(irep = 1:250) %dorng% sample_coupled_chains(mixkernel, coupled_mixkernel, init_chains, tuning, m = 1, lag = 1)

## plot histogram of meeting times
meeting_times <- sapply(coupled_chains_tilmeet, function(l) l$meetingtime)
qplot(x = meeting_times, geom = "histogram") + theme_minimal() + xlab("meeting time") +
  theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20))




