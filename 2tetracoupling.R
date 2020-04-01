rm(list=ls())

library(ggplot2)
theme_set(theme_void())
library(gridExtra)

colors <- c(rgb(1, 0.1, 0.3), rgb(0.3, 0.1, 1))

## various couplings

# normal means
mu1 <- -1
mu2 <- 2
# std deviation
sigma <- .7
# number of samples
nsamples <- 1e4
# reflection-maximal coupling
reflmax_samples <- matrix(nrow=nsamples, ncol=2)
# draw x components first
xdot <- rnorm(nsamples)
# this follows the notation of Bou Rabee et al, 2018, roughly
z <- (mu1 - mu2) / sigma
normz <- sqrt(sum(z^2))
e <- z / normz
utilde <- runif(nsamples, 0, 1)
accepts <- (log(utilde) < (dnorm(xdot + z, 0, 1, log=TRUE)
                          - dnorm(xdot, log=TRUE)))
ydot <- rep(0, nsamples)
ydot[accepts] <- (xdot)[accepts] + z
ydot[!accepts] <-  xdot[!accepts] - 2 * (e * xdot[!accepts]) * e
reflmax_samples[, 1] <- mu1 + sigma * xdot
reflmax_samples[, 2] <- mu2 + sigma * ydot

df <- data.frame(coupling=rep("reflection-maximal", nsamples),
                 x=reflmax_samples[, 1],
                 y=reflmax_samples[, 2])

gscatter <- qplot(x=reflmax_samples[, 1], y=reflmax_samples[, 2],
                  geom="blank") +
            geom_point(alpha=0.15) +
            geom_abline(slope=1, intercept=0, linetype=2)
gmargx <- qplot(x=reflmax_samples[, 1], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[1]) +
            stat_function(fun=function(x) dnorm(x, mu1, sigma))
gmargy <- qplot(x=reflmax_samples[, 2], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[2]) +
            stat_function(fun=function(x) dnorm(x, mu2, sigma)) +
            coord_flip() +
            scale_y_reverse()
empty <- ggplot()

g <- gridExtra::grid.arrange(empty, gmargx, gmargy, gscatter,
                            ncol=2, nrow=2,
                            widths=c(1, 4), heights=c(1, 4))
# ggsave(filename="variouscouplings.1.pdf", plot=g)

# reflection coupling
refl_samples <- matrix(0, nrow=nsamples, ncol=2)
refl_samples[, 1] <- reflmax_samples[, 1]
refl_samples[, 2] <- mu2 - (refl_samples[, 1] - mu1)

gscatter <- qplot(x=refl_samples[, 1], y=refl_samples[, 2], geom="blank") +
            geom_point(alpha=0.15) +
            geom_abline(slope=1, intercept=0, linetype=2)
gmargx <- qplot(x=refl_samples[, 1], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[1]) +
            stat_function(fun=function(x) dnorm(x, mu1, sigma))
gmargy <- qplot(x=refl_samples[, 2], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[2]) +
            stat_function(fun=function(x) dnorm(x, mu2, sigma)) +
            coord_flip()
empty <- ggplot()

g <- gridExtra::grid.arrange(gmargx, empty, gscatter, gmargy,
                            ncol=2, nrow=2,
                            widths=c(4, 1), heights=c(1, 4))
# ggsave(filename="variouscouplings.2.pdf", plot=g)

# optimal transport coupling
transport_samples <- matrix(0, nrow=nsamples, ncol=2)
transport_samples[, 1] <- reflmax_samples[, 1]
transport_samples[, 2] <- mu2 - mu1 + transport_samples[, 1]

gscatter <- qplot(x=transport_samples[, 1], y=transport_samples[, 2],
                 geom="blank") +
            geom_point(alpha=0.15) +
            geom_abline(slope=1, intercept=0, linetype=2)
gmargx <- qplot(x=transport_samples[, 1], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[1]) +
            stat_function(fun=function(x) dnorm(x, mu1, sigma)) +
            scale_y_reverse()
gmargy <- qplot(x=transport_samples[, 2], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[2]) +
            stat_function(fun=function(x) dnorm(x, mu2, sigma)) +
            coord_flip() +
            scale_y_reverse()
empty <- ggplot()

g <- gridExtra::grid.arrange(gmargy, gscatter, empty, gmargx,
                          ncol=2, nrow=2,
                          widths=c(1,4), heights=c(4, 1))
# ggsave(filename="variouscouplings.3.pdf", plot=g)

# max coupling
max_samples <- matrix(0, nrow=nsamples, ncol=2)
max_samples[, 1] <- reflmax_samples[, 1]
dp <- function(x) dnorm(x, mean=mu1, sd=sigma, log=TRUE)
dq <- function(x) dnorm(x, mean=mu2, sd=sigma, log=TRUE)
rq <- function(n) rnorm(n, mean=mu2, sd=sigma)
for (isample in 1:nsamples){
	x <- max_samples[isample, 1]
	if (dp(x) + log(runif(1)) < dq(x)){
		max_samples[isample, 2] <- x
	} else {
		reject <- TRUE
		y <- NA
		while (reject){
			y <- rq(1)
			reject <- (dq(y) + log(runif(1)) < dp(y))
		}
		max_samples[isample, 2] <- y
	}
}

gscatter <- qplot(x=max_samples[, 1], y=max_samples[, 2], geom="blank") +
            geom_point(alpha=0.15) +
            geom_abline(slope=1, intercept=0, linetype=2)
gmargx <- qplot(x=max_samples[, 1], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[1]) +
            stat_function(fun=function(x) dnorm(x, mu1, sigma)) +
            scale_y_reverse()
gmargy <- qplot(x=max_samples[, 2], geom="blank") +
            geom_histogram(aes(y=..density..), fill=colors[2]) +
            stat_function(fun=function(x) dnorm(x, mu2, sigma)) +
            coord_flip()
empty <- ggplot()

g <- gridExtra::grid.arrange(gscatter, gmargy, gmargx, empty,
                            ncol=2, nrow=2,
                            widths=c(4, 1), heights=c(4, 1))
# ggsave(filename="variouscouplings.4.pdf", plot=g)
