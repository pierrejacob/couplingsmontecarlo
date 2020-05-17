##
## transition matrix
rm(list = ls())
# source("couplingsmontecarlo/commongrounds.R")
library(doParallel)
library(doRNG)
registerDoParallel(cores = 10)
library(ggplot2)
library(ggridges)
library(reshape2)
library(dplyr)
set.seed(6)
set_theme_chapte5 <- function(){
  library(ggplot2)
  library(gridExtra)
  theme_set(theme_void())
  colors <- c(rgb(0.1,0.7,0.4), rgb(0.1, 0.35, 0.35))
  return(list(colors = colors))
}
graphsettings <- set_theme_chapte5()

nstates <- 5

transmatrix <- matrix(NA, nrow = nstates, ncol = nstates)
for (istate in 1:nstates){
  transmatrix[istate,] <- rexp(nstates)
  transmatrix[istate,] <- transmatrix[istate,] / sum(transmatrix[istate,])
}
## find stationary distribution
invar_vec <- eigen(t(transmatrix))$vectors[,1]
invar_vec <- invar_vec / sum(invar_vec)
invar_vec %*% transmatrix

## random function representation
transmap <- function(x,u) 1 + findInterval(u, cumsum(transmatrix[x,]))

# ## represent function with arrows
# from_ <- 1:nstates
# u_ <- runif(1)
# to_ <- sapply(from_, function(x) transmap(x, u_))
# data.frame(time = , from_ = from_, to_ = to_)

## plot with arrows
set.seed(1)
times_ <- (-4):0
nsteps <- length(times_)
us <- runif(nsteps)
df_from_to <- data.frame()
for (step in 1:(nsteps-1)){
  to_ <- sapply(1:nstates, function(x) transmap(x, us[step]))
  df_from_to <- rbind(df_from_to, data.frame(time = times_[step], from = 1:nstates, to = to_))
}
head(df_from_to)

gmaps <- ggplot(df_from_to, aes(x = 0, xend = 1, y = from, yend = to)) + geom_segment(arrow = arrow(angle = 5, type = 'closed', length = unit(.5, "cm")))
gmaps <- gmaps + geom_point(size = 5, shape = 1) + xlab('') + ylab('')
gmaps <- gmaps + coord_flip() + facet_wrap(~ time, nrow = 1) + theme(strip.background = element_blank(), strip.text.x = element_blank())
gmaps
ggsave(filename = "../5maps.pdf", plot = gmaps, width = 10, height = 5)
#
g <- ggplot(df_from_to, aes(x = time, xend = time + 1, y = from, yend = to))
g <- g + geom_point(size = 5, shape = 1) + xlab('') + ylab('')
g <- g + geom_point(data = data.frame(x = rep(0, nstates), y = 1:nstates),
                    aes(x = x, y = y, xend = NULL, yend = NULL), size = 5, shape = 1)
g <- g + geom_segment(arrow = arrow(angle = 5, type = 'closed', length = unit(.5, "cm")))
g <- g + scale_x_continuous(breaks = times_) + scale_y_continuous(breaks = 1:nstates)
g <- g + theme_test() + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
g <- g + theme(panel.grid.major.x = element_line(colour = 'grey'))
g

ggsave(filename = "../5cftp.pdf", plot = g, width = 10, height = 5)
#
# ##
# ## CFTP
# ## with time decreasing by one
# CFTP <- function(transmap){
#   coalescence <- FALSE
#   time_ <- 0
#   u <- c()
#   # store the composition of a certain number of maps
#   composition <- 1:nstates
#   while (!(coalescence)){
#     time_ <- time_ - 1
#     u <- c(runif(1), u)
#     ## starting from anywhere
#     xs <- 1:nstates
#     ## apply new random function
#     xs <- sapply(xs, function(x) transmap(x, u[1]))
#     ## apply previous composition
#     xs <- composition[xs]
#     ## update composition
#     composition <- xs
#     # check coalescence
#     coalescence <- (length(unique(xs))==1)
#   }
#   return(list(time = time_, x = unique(xs)))
# }
#
# nrep <- 1e4
# xs <- lapply(1:nrep, function(i) CFTP(transmap))
# table(sapply(xs, function(z) z$x))/nrep
# invar_vec
# hist(sapply(xs, function(z) z$time))
#
# ## what is the point of this doubling mechanism?
#
# CFTP_doubling <- function(transmap){
#   coalescence <- FALSE
#   time_ <- -1
#   u <- c()
#   while (!(coalescence)){
#     if (time_ == -1){
#       u <- runif(n = 1)
#     } else {
#       u <- c(runif(n = abs(time_)/2), u)
#     }
#     # check coalescence
#     outputstates <- sapply(1:nstates, function(x){
#       for (i in 1:length(u)){
#         x <- transmap(x, u[i])
#       }
#       return(x)
#     })
#     coalescence <- (length(unique(outputstates))==1)
#     time_ <- 2 * time_
#   }
#   return(list(time = time_, x = unique(outputstates)))
# }
#
# nrep <- 1e4
# xs <- lapply(1:nrep, function(i) CFTP_doubling(transmap))
# table(sapply(xs, function(z) z$x))/nrep
# invar_vec
# hist(sapply(xs, function(z) z$time))
#
# ## we could check that state at forward coalescence not distributed according
# ## to pi
