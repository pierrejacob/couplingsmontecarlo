
rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter4()
library(reshape2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

tv_upper_bound <- function(meetingtimes, lag, t){
    return(mean(pmax(0, ceiling((meetingtimes-lag-t) / lag))))
}
rmax_coupling <- function(rp, dp, rq, dq){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)) {
        return(list(xy = c(x, x), identical = TRUE))
    }
    else {
        reject <- TRUE
        y <- NA
        while (reject) {
            y <- rq(1)
            reject <- (dq(y) + log(runif(1)) < dp(y))
        }
        return(list(xy = c(x, y), identical = FALSE))
    }
}

## Compute log-density of inverse Gamma
## at x > 0 and with given parameters alpha, beta, given by
##  alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
dinversegamma <- function(x, alpha, beta){
    return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

## Sample from inverse Gamma with parametrization as above
rinversegamma <- function(n, alpha, beta){
    return(1/rgamma(n = n, shape = alpha, rate = beta))
}

## Sample from maximally coupled inverse Gamma variables
rinversegamma_coupled <- function(alpha1, alpha2, beta1, beta2){
    return(rmax_coupling(function(n) rinversegamma(n, alpha1, beta1),
                         function(x) dinversegamma(x, alpha1, beta1),
                         function(n) rinversegamma(n, alpha2, beta2),
                         function(x) dinversegamma(x, alpha2, beta2)))
}

## Sample from maximally coupled Normal variables (a reflection-maximal coupling might be better but hey).
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
    return(rmax_coupling(function(n) rnorm(n, mu1, sigma1),
                         function(x) dnorm(x, mu1, sigma1, log = TRUE),
                         function(n) rnorm(n, mu2, sigma2),
                         function(x) dnorm(x, mu2, sigma2, log = TRUE)))
}
## the Dyestuff data set is available in the package lme4
library(lme4)
# ?Dyestuff
mu_0 <- 0; sigma_0 <- 10^6; sigma_02 <- sigma_0^2
a1 <- 0.5; b1 <- 1
a2 <- 0; b2 <- 0
K <- 6; J <- 5
Ydataframe <- Dyestuff %>% group_by(Batch) %>% mutate(individual = 1:n()) %>% ungroup()
# let's look at the data
Ydataframe %>% group_by(Batch) %>% summarise(mean_per_batch = mean(Yield), var_per_batch = var(Yield))
# manipulate the data to ease forthcoming calculations
Ywide <- tidyr::spread(Ydataframe, Batch, Yield) %>% select(-individual)
Ywide <- t(as.matrix(Ywide)) # data as matrix, with each row corresponding to a batch
Ymean_per_batch <- rowMeans(Ywide)
Ymean <- mean(Ywide)
# note: here we're doing some really advanced stuff: removing a vector of length K to a matrix with K rows
# does a row-wise subtraction with the elements of the vector
# e.g. Ywide - c(1,1,1,2,2,2) removes 1 to the first 3 rows, and removes 2 to the next 3 rows

## Gibbs sampler, as in the CR paper
single_kernel <- function(state){
    # extract parameters from vector state
    sigma_theta2 <- state[1];  sigma_e2 <- state[2];  mu <- state[3];  theta <- state[4:(4+K-1)]
    # update of sigma_theta2
    sigma_theta2 <- rinversegamma(1, a1 + 0.5 * K, b1 + 0.5 * sum((theta - mu)^2))
    # update of sigma_22
    sigma_e2 <- rinversegamma(1, a2 + 0.5 * K * J, b2 + 0.5 * sum((Ywide - theta)^2))
    # update of mu
    mean_mu <- (sigma_theta2 * mu_0 + sigma_02 * sum(theta)) / (sigma_theta2 + K * sigma_02)
    var_mu <- (sigma_theta2*sigma_02)/(sigma_theta2 + K * sigma_02)
    mu <- rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
    # update of each theta
    mean_theta <- (J * sigma_theta2 * Ymean_per_batch + sigma_e2 * mu) / (J * sigma_theta2 + sigma_e2)
    var_theta <- (sigma_theta2 * sigma_e2) / (J * sigma_theta2 + sigma_e2)
    theta <- rnorm(K, mean = mean_theta, sd = sqrt(var_theta))
    return(c(sigma_theta2, sigma_e2, mu, theta))
}


## Coupled Gibbs sampler, with maximally coupled updates
coupled_kernel <- function(state1, state2){
    sigma_theta21 <- state1[1];  sigma_e21 <- state1[2];  mu1 <- state1[3];  theta1 <- state1[4:(4+K-1)]
    sigma_theta22 <- state2[1];  sigma_e22 <- state2[2];  mu2 <- state2[3];  theta2 <- state2[4:(4+K-1)]
    whichidentical <- rep(FALSE, length(state1))
    identical <- TRUE
    # update of sigma_theta2
    sigma_theta2_ <- rinversegamma_coupled(a1 + 0.5 * K, a1 + 0.5 * K,
                                           b1 + 0.5 * sum((theta1 - mu1)^2),
                                           b1 + 0.5 * sum((theta2 - mu2)^2))
    sigma_theta21 <- sigma_theta2_$xy[1]; sigma_theta22 <- sigma_theta2_$xy[2]
    whichidentical[1] <- sigma_theta2_$identical
    # update of sigma_22
    sigma_e2_ <- rinversegamma_coupled(a2 + 0.5 * K * J, a2 + 0.5 * K * J,
                                       b2 + 0.5 * sum((Ywide - theta1)^2),
                                       b2 + 0.5 * sum((Ywide - theta2)^2))
    sigma_e21 <- sigma_e2_$xy[1]; sigma_e22 <- sigma_e2_$xy[2]
    whichidentical[2] <- sigma_e2_$identical
    # update of mu
    mean_mu1 <- (sigma_theta21 * mu_0 + sigma_02 * sum(theta1)) / (sigma_theta21 + K * sigma_02)
    var_mu1  <- (sigma_theta21*sigma_02)/(sigma_theta21 + K * sigma_02)
    mean_mu2 <- (sigma_theta22 * mu_0 + sigma_02 * sum(theta2)) / (sigma_theta22 + K * sigma_02)
    var_mu2  <- (sigma_theta22*sigma_02)/(sigma_theta22 + K * sigma_02)
    mu_ <- rnorm_max_coupling(mean_mu1, mean_mu2, sqrt(var_mu1), sqrt(var_mu2))
    mu1 <- mu_$xy[1]; mu2 <- mu_$xy[2]
    whichidentical[3] <- mu_$identical
    # update of each theta
    mean_theta1 <- (J * sigma_theta21 * Ymean_per_batch + sigma_e21 * mu1) / (J * sigma_theta21 + sigma_e21)
    var_theta1 <- (sigma_theta21 * sigma_e21) / (J * sigma_theta21 + sigma_e21)
    mean_theta2 <- (J * sigma_theta22 * Ymean_per_batch + sigma_e22 * mu2) / (J * sigma_theta22 + sigma_e22)
    var_theta2 <- (sigma_theta22 * sigma_e22) / (J * sigma_theta22 + sigma_e22)
    for (k in 1:K){
        theta_k <- rnorm_max_coupling(mean_theta1[k], mean_theta2[k], sqrt(var_theta1), sqrt(var_theta2))
        theta1[k] <- theta_k$xy[1]; theta2[k] <- theta_k$xy[2]
        whichidentical[3+k] <- theta_k$identical
    }
    return(list(state1 = c(sigma_theta21, sigma_e21, mu1, theta1),
                state2 = c(sigma_theta22, sigma_e22, mu2, theta2),
                identical = all(whichidentical), whichidentical = whichidentical))
}
v1 <- mean((Ywide - Ymean_per_batch)^2)
v2 <- mean((Ymean_per_batch - Ymean)^2)
rinit <- function() c(1, 1, Ymean, (J * v1 * Ymean_per_batch + v2 * Ymean) / (J * v1 + v2))


sample_coupled_chains <- function (single_kernel, coupled_kernel, rinit, m = 1, lag = 1,
                                   max_iterations = Inf, preallocate = 10)
{
    starttime <- Sys.time()
    state1 <- rinit()
    state2 <- rinit()
    dimstate <- length(state1)
    nrowsamples1 <- m + preallocate + lag
    samples1 <- matrix(nrow = nrowsamples1, ncol = dimstate)
    samples2 <- matrix(nrow = nrowsamples1 - lag, ncol = dimstate)
    idcomponents <- matrix(nrow = nrowsamples1 - lag, ncol = dimstate)
    samples1[1, ] <- state1
    samples2[1, ] <- state2
    idcomponents[1,] <- rep(FALSE, dimstate)
    time <- 0
    for (t in 1:lag) {
        time <- time + 1
        state1 <- single_kernel(state1)
        samples1[time + 1, ] <- state1
    }
    meetingtime <- Inf
    while ((time < max(meetingtime, m)) && (time < max_iterations)) {
        time <- time + 1
        if (is.finite(meetingtime)) {
            state1 <- single_kernel(state1)
            state2 <- state1
            idc_ <- rep(TRUE, dimstate)
        } else {
            res_coupled_kernel <- coupled_kernel(state1, state2)
            state1 <- res_coupled_kernel$state1
            state2 <- res_coupled_kernel$state2
            if (res_coupled_kernel$identical){
                meetingtime <- time
            }
            idc_ <- res_coupled_kernel$whichidentical
        }
        if ((time + 1) > nrowsamples1) {
            new_rows <- nrowsamples1
            nrowsamples1 <- nrowsamples1 + new_rows
            samples1 <- rbind(samples1, matrix(NA, nrow = new_rows,
                                               ncol = dimstate))
            samples2 <- rbind(samples2, matrix(NA, nrow = new_rows,
                                               ncol = dimstate))
            idcomponents <- rbind(idcomponents, matrix(NA, nrow = new_rows,
                                               ncol = dimstate))
        }
        samples1[time + 1, ] <- state1
        samples2[time - lag + 1, ] <- state2
        idcomponents[time - lag + 1,] <- idc_
    }
    samples1 <- samples1[1:(time + 1), , drop = F]
    samples2 <- samples2[1:(time - lag + 1), , drop = F]
    idcomponents <- idcomponents[1:(time - lag + 1), , drop = F]
    cost <- lag + 2 * (meetingtime - lag) + max(0, time - meetingtime)
    currenttime <- Sys.time()
    elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) -
                                                         lubridate::ymd_hms(starttime)), "seconds")
    return(list(samples1 = samples1, samples2 = samples2, meetingtime = meetingtime,
                iteration = time, elapsedtime = elapsedtime, cost = cost, idcomponents = idcomponents))
}

lag <- 300

##


## Now, how many steps do we need to be close to stationarity in TV??
NREP <- 5e2
coupledchains <- foreach(irep = 1:NREP) %dorng% { sample_coupled_chains(single_kernel, coupled_kernel, rinit, lag = lag)}
meeting_times <- sapply(coupledchains, function(c) c$meetingtime)
niterations <- 5e2
upperbounds <- sapply(1:niterations, function(t) tv_upper_bound(meeting_times, lag, t))
# gridExtra::grid.arrange(ghist, g_tvbounds + theme_bw(), nrow = 1)

hist_noplot <- hist(meeting_times - lag, plot = F, nclass = 50)
xgrid <- c(min(hist_noplot$breaks), hist_noplot$mids, max(hist_noplot$breaks))
densitygrid <- c(0, hist_noplot$density, 0)


g_tvbounds <- qplot(x = 1:niterations, y = upperbounds, geom = "line")
g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
g_tvbounds <- g_tvbounds + scale_y_continuous(breaks = c(0,1), limits = c(0,1.1))
g_tvbounds <- g_tvbounds + geom_ribbon(data=data.frame(x = xgrid,
                                                       ymin = rep(0, length(xgrid)),
                                                       y = densitygrid/max(densitygrid)),
                                       aes(x= x, ymin = ymin, ymax = y, y=NULL), alpha = .3, fill = graphsettings$colors[2]) + geom_line()
g_tvbounds <- g_tvbounds + theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
                                 axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20, angle = 90))

pdf("../gibbs.tvbounds.pdf", width = 10, height = 5)
print(g_tvbounds)
dev.off()


index <- sample(x = which(meetingtimes-lag>100), size = 1)
cc <- coupledchains[[index]]
idcomponents.df <- lapply(as.list(1:ncol(cc$idcomponents)), function(l){
    data.frame(id = cc$idcomponents[,l], time = 1:nrow(cc$idcomponents), component = l)
}) %>% bind_rows()
gcomponents <- ggplot(idcomponents.df, aes(x = time, y = component, fill = factor(id))) + geom_tile() +
    scale_fill_manual(values=c("white", graphsettings$colors[2])) +
    theme(legend.position="none")
library(latex2exp)
label.df <- data.frame(x = -3, y = 1:(3+K), label = c("$\\sigma_\\theta$", "$\\sigma_e$", "$\\mu$", paste0("$\\theta_", 1:K, "$")))
gcomponents <- gcomponents + annotate("text", x = label.df$x, y = label.df$y,
                       label = lapply(label.df$label, function(x){TeX(x, output = "character")}), parse = T, size = 5)
gcomponents

ggsave(filename = "../4gibbscomponents.pdf", plot = gcomponents, width = 10, height = 5)
