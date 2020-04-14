rm(list = ls())
source("commongrounds.R")
library(doParallel)
library(doRNG)
registerDoParallel(cores = 10)

# 
horizon <- 200
arcoef <- 0.95
sigma <- 1
xchain <- rep(0, horizon+1)
ychain <- rep(0, horizon)

xchain[1] <- rnorm(1) 
xchain[2] <- arcoef * xchain[1] + rnorm(1, sd = sigma)
ychain[1] <- rnorm(1)
for (t in 2:horizon){
  increment <- rnorm(1, mean = 0, sd = sigma)
  xchain[t+1] <- arcoef * xchain[t] + increment
  ychain[t] <- arcoef * ychain[t-1] + increment
}
par(mfrow = c(2,1))
matplot(cbind(xchain, ychain), type = "l", ylab = "chains", xlab = "time")
matplot(abs(xchain[2:horizon] - ychain[1:(horizon-1)]), type = "l", ylab = "chains", xlab = "time")

