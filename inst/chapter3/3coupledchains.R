rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
set_theme_chapter3 <- function(){
  library(ggplot2)
  library(gridExtra)
  theme_set(theme_void())
  colors <- c(rgb(0.8,0.5,0.2), rgb(0.2, 0.6, 0.9))
  return(list(colors = colors))
}
graphsettings <- set_theme_chapter3()
# library(ggridges)
# library(reshape2)
# library(dplyr)



# logtarget <- function(x) log(0.5 * dnorm(x, mean = 0, sd = 1) + 0.5 * dnorm(x, mean = -4, sd = 1))
# logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
logtarget <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 4), sd = c(1.2,2.1), log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
curve(sapply(x, function(xp) exp(logtarget(xp))), from = -5, to = 8)
sd_proposal <- 1
rinit <- function() rnorm(1, 5, 1)

MH_kernel <- function(chain_state){
  proposal_value <- rnorm(1, chain_state, sd_proposal)
  proposal_pdf <- logtarget(proposal_value)
  current_pdf <- logtarget(chain_state)
  if (log(runif(1)) < (proposal_pdf - current_pdf)){
    return(proposal_value)
  } else {
    return(chain_state)
  }
}

coupled_MH_kernel <- function(chain_state1, chain_state2){
  maxnormal <- getmaxcoupling(function(n) rnorm(n, chain_state1, sd_proposal),
                              function(x) dnorm(x, chain_state1, sd_proposal, log = TRUE),
                              function(n) rnorm(n, chain_state2, sd_proposal),
                              function(x) dnorm(x, chain_state2, sd_proposal, log = TRUE))
  proposal_value <- maxnormal()
  proposal_pdf1 <- logtarget(proposal_value$xy[1])
  proposal_pdf2 <- logtarget(proposal_value$xy[2])
  current_pdf1 <- logtarget(chain_state1)
  current_pdf2 <- logtarget(chain_state2)
  logu <- log(runif(1))
  if (is.finite(proposal_pdf1)){
    if (logu < (proposal_pdf1 - current_pdf1)){
      chain_state1 <- proposal_value$xy[1]
    }
  }
  if (is.finite(proposal_pdf2)){
    if(logu < (proposal_pdf2 - current_pdf2)){
      chain_state2 <- proposal_value$xy[2]
    }
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
}

set.seed(1)
nmcmc <- 100
xchain <- rep(6, nmcmc)
ychain <- rep(-2, nmcmc)

for (imcmc in 2:nmcmc){
  res_ <- coupled_MH_kernel(xchain[imcmc-1], ychain[imcmc-1])
  xchain[imcmc] <- res_$chain_state1
  ychain[imcmc] <- res_$chain_state2
}

g <- qplot(x = 1:nmcmc, y = xchain, geom = "blank") + geom_line(colour = graphsettings$colors[1]) +
  geom_line(aes(y = ychain), colour = graphsettings$colors[2])
g
ggsave(filename = "../coupledmcmc.path.pdf", plot = g, width = 10, height = 7)

