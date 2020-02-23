rm(list = ls())
source("commongrounds.R")

p <- function(x) dbeta(x, 5, 5)
q <- function(x) dbeta(x, 2, 1.5)

rmc <- getmaxcoupling(rp = function(n) rbeta(n, 5, 5), dp = function(x) dbeta(x, 5, 5, log = TRUE),
                      rq = function(n) rbeta(n, 2, 1.5), dq = function(x) dbeta(x, 2, 1.5, log = TRUE))

res <- sapply(1:1e4, function(i) rmc()$xy)

gscatter <- qplot(x = res[1,], y = res[2,], geom = "blank") + geom_point(alpha = 0.5) 
gmargx <- qplot(x = res[1,], geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[1]) + stat_function(fun = p)
gmargy <- qplot(x = res[2,], geom = "blank") + geom_histogram(aes(y=..density..), fill = colors[2]) + stat_function(fun = q) + coord_flip()
empty <- ggplot()

g <- gridExtra::grid.arrange(gmargx, empty, gscatter, gmargy, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
g

# ggsave(filename = "tvdistance2.pdf", plot = g3)
