rm(list = ls())
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter3()
set.seed(1)

n <- 25
x <- matrix(0.1, nrow = n, ncol = n)
x[5:15,5:15] <- 0.5
x[3:10,3:10] <- 0.9
x[15:n,10:n] <- 0.75

image(x, col = gray.colors(10))
df <- expand.grid(1:n,1:n)
df$z <- as.numeric(x)
gx <- ggplot(df, aes(x=Var1, y=Var2, fill=z)) +
    geom_tile() + scale_fill_gradient(low = "black", high = "white") + theme(legend.position = "none")
gx
# ggsave(filename = "../imagerestoration1.pdf", gx)


# library(truncnorm)
# y <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = x, sd = 1e-5)
y <- matrix(0, n, n)
for (i in 1:n){
    for (j in 1:n){
        y[i,j] <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = x[i,j], sd = 2e-1)
        # y[i,j] <- TruncatedNormal::rtnorm(n = 1, mu = x[i,j], sd = 1e-1, lb = 0, ub = 1)

    }
}
df$z <- as.numeric(y)
gy <- ggplot(df, aes(x=Var1, y=Var2, fill=z)) +
    geom_tile() + scale_fill_gradient(low = "black", high = "white") + theme(legend.position = "none")
gy
ggsave(filename = "../imagerestoration2.pdf", gy)

##
rinit <- function(){
    return(matrix(runif(n*n), nrow = n, ncol = n))
}
##
state <- rinit()
gamma2 <- (4)^2
sigma <- 0.2
precision <- 1/(sigma^2)
single_kernel <- function(state){
    ni <- 4
    for (i in 1:n){
        for (j in 1:n){
            cond_var <- 1/(precision + gamma2 * ni)
            sum_neighbors <- state[i, ifelse(j < n, j + 1, 1)] + state[i, ifelse(j > 1, j - 1, n)] +
                state[ifelse(i > 1, i - 1, n), j] + state[ifelse(i < n, i + 1, 1), j]
            cond_mean <- cond_var * (precision * y[i,j] + gamma2 * sum_neighbors)
            state[i,j] <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = cond_mean, sd = sqrt(cond_var))
        }
    }
    return(state)
}

nmcmc <- xe3
burnin <- 1e3
state_history <- array(NA, dim = c(n, n, nmcmc))
state_sum <- matrix(0, n, n)
for (imcmc in 1:nmcmc){
    state <- single_kernel(state)
    state_history[,,imcmc] <- state
    if (imcmc > burnin){
        state_sum <- state_sum + state
    }
}

df$z <- as.numeric(state_sum/(nmcmc-burnin))
grestore <- ggplot(df, aes(x=Var1, y=Var2, fill=z)) +
    geom_tile() + scale_fill_gradient(low = "black", high = "white") + theme(legend.position = "none")
grestore
ggsave(filename = "../imagerestoration3.pdf", grestore)

