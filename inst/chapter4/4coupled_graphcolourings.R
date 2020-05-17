## implementation of coupled Gibbs sampler on graph colourings
rm(list=ls())
set.seed(1)
library(couplingsmontecarlo)
graphsettings <- set_theme_chapter4()
library(ggridges)
library(reshape2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
# load igraph package to deal with graphs
library(igraph)
# function to plot graph
custom_plot_graph <- function(g) plot(g, layout=layout_with_kk, vertex.label = NA)
custom_plot_2graphs <- function(g1, g2){
    bivalues <- lapply(1:nvertices, function(index) c(1,1))
    bicols <- lapply(1:nvertices, function(index) c(V(g1)$color[index], V(g2)$color[index]))
    plot(g1, layout = layout_with_kk, vertex.shape="pie", vertex.pie=bivalues, vertex.pie.color=bicols,  vertex.label=NA)
}

## create graph
graph_name <- "Icosahedral"
g <- igraph::make_graph(graph_name)
## number of vertices in the graph
nvertices <- length(V(g))
## maximal degree of graph
maxdegree <- max(igraph::degree(g))
## number of colours
q <- 2 * maxdegree + 2
## all possible colours
all_colours <- colorRampPalette(RColorBrewer::brewer.pal(11, name = "PuOr"))(q)

## greedy colouring to initialize chains
rinit <- function(g){
    V(g)$color <- "black"
    V(g)$color[1] <- all_colours[1]
    already_colored <- c(1)
    for (ivertex in 2:nvertices){
        ## get neighbors
        n_i <- neighbors(g, ivertex)
        ## get colors not on neighbouring vertices
        legal_colours <- setdiff(all_colours, V(g)$color[n_i])
        ## color considered vertex with first legal colour
        V(g)$color[ivertex] <- legal_colours[1]
    }
    return(g)
}

## Markov chain, Gibbs sampler on graph colourings
single_kernel <- function(g){
    ## choose a vertex at random
    ivertex <- sample(1:nvertices, 1)
    ## get neighbours
    n_i <- neighbors(g, ivertex)
    ## get colors not in neighbors
    legal_colours <- setdiff(all_colours, V(g)$color[n_i])
    ## randomly sample one such colour
    V(g)$color[ivertex] <- sample(x = legal_colours, size = 1)
    return(g)
}


## coupled Gibbs sampler kernel
##
coupled_kernel <- function(g1, g2){
    # choose a vertex at random
    ivertex <- sample(1:nvertices, 1)
    V(g1)$color[ivertex]
    V(g2)$color[ivertex]
    # find neighbors
    n_i <- neighbors(g1, ivertex)
    # get colors not in neighbors
    legal_colours1 <- setdiff(all_colours, V(g1)$color[n_i])
    legal_colours2 <- setdiff(all_colours, V(g2)$color[n_i])
    ## sample max coupling of unif distributions over legal colours
    p1 <- rep(0, q)
    p2 <- rep(0, q)
    p1[match(legal_colours1, all_colours)] <- 1/length(legal_colours1)
    p2[match(legal_colours2, all_colours)] <- 1/length(legal_colours2)
    # decompose into common part
    common_part <- pmin(p1, p2)
    c <- sum(common_part)
    # and residual parts
    r1 <- (p1-common_part)/(1-c)
    r2 <- (p2-common_part)/(1-c)
    # sample pair of indices
    if (runif(1) < c){
        index1 <- sample(x = 1:q, size = 1, prob = common_part/c)
        index2 <- index1
    } else {
        index1 <- sample(x = 1:q, size = 1, prob = r1)
        index2 <- sample(x = 1:q, size = 1, prob = r2)
    }
    ## colour nodes
    V(g1)$color[ivertex] <- all_colours[index1]
    V(g2)$color[ivertex] <- all_colours[index2]
    ## return two graphs
    return(list(g1 = g1, g2 = g2))
}


## initialize chains
g1 <- rinit(g)
g2 <- rinit(g)
## advance one chain for 'lag' steps
lag <- 1e3
for (iter in 1:lag){
    g1 <- single_kernel(g1)
}
## plot graphs with two colours per node
#pdf("../twographs.pdf")
custom_plot_2graphs(g1, g2)
#dev.off()

## construction of two chains with a lag L
## using single_kernel and coupled_kernel
sample_meetingtime <- function(single_kernel, coupled_kernel, rinit, lag = 1, max_iterations = Inf){
    starttime <- Sys.time()
    # initialize two chains
    state1 <- rinit(g); state2 <- rinit(g)
    # move first chain for 'lag' iterations
    time <- 0
    for (t in 1:lag){
        time <- time + 1
        state1 <- single_kernel(state1)
    }
    # move two chains until meeting (or until max_iterations)
    meetingtime <- Inf
    # two chains could already be identical by chance
    if (all(V(state1)$color == V(state2)$color)) meetingtime <- lag
    while (is.infinite(meetingtime) && (time < max_iterations)){
        time <- time + 1
        # use coupled kernel
        coupledstates <- coupled_kernel(state1, state2)
        state1 <- coupledstates$g1
        state2 <- coupledstates$g2
        # check if meeting has occurred
        if (all(V(state1)$color == V(state2)$color)) meetingtime <- time
    }
    currentime <- Sys.time()
    elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
    return(list(meetingtime = meetingtime, elapsedtime = elapsedtime))
}

## generate a number of meeting times, for a certain lag
nrep <- 1e3
lag <- 200

sample_meetingtime(single_kernel, coupled_kernel, rinit, lag)

meetings <- foreach(irep = 1:nrep) %dorng% sample_meetingtime(single_kernel, coupled_kernel, rinit, lag, max_iterations = 1e5)
meeting_times <- sapply(meetings, function(x) x$meetingtime)
## plot meeting times
ghist <- qplot(x = meeting_times - lag, geom = "blank") + geom_histogram(aes(y=..density..))
ghist <- ghist + xlab("meeting time - lag")
ghist <- ghist + theme_minimal()
ghist
## compute TV upper bounds
tv_upper_bound_estimates <- function(meeting_times, L, t){
    return(mean(pmax(0,ceiling((meeting_times-L-t)/L))))
}
niter <- (floor(1.1*max(meeting_times)-lag))
upperbounds <- sapply(1:niter, function(t) tv_upper_bound_estimates(unlist(meeting_times), lag, t))
g_tvbounds <- qplot(x = 1:niter, y = upperbounds, geom = "line")
g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
# g_tvbounds <- g_tvbounds + labs(title = paste0("uniform colorings of ", graph_name, " graph"))
g_tvbounds <- g_tvbounds + scale_y_continuous(breaks = (1:10)/10, limits = c(0,1.1))
g_tvbounds <- g_tvbounds + theme_minimal()
gridExtra::grid.arrange(ghist, g_tvbounds, nrow = 1)

# pdf("../graphcolourings.meetingtimes.pdf")
# ghist
# dev.off()

# pdf("../graphcolourings.tvbounds.pdf")
# g_tvbounds
# dev.off()

hist_noplot <- hist(meeting_times - lag, plot = F, nclass = 30)
xgrid <- c(min(hist_noplot$breaks), hist_noplot$mids, max(hist_noplot$breaks))
densitygrid <- c(0, hist_noplot$density, 0)
g_tvbounds <- qplot(x = 1:niter, y = upperbounds, geom = "line")
g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
g_tvbounds <- g_tvbounds + scale_y_continuous(breaks = c(0,1), limits = c(0,1.1))
g_tvbounds <- g_tvbounds + geom_ribbon(data=data.frame(x = xgrid,
                                         ymin = rep(0, length(xgrid)),
                                         y = densitygrid/max(densitygrid)),
                         aes(x= x, ymin = ymin, ymax = y, y=NULL), alpha = .3, fill = graphsettings$colors[2]) + geom_line()
g_tvbounds <- g_tvbounds + theme(axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
                                 axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20, angle = 90))

pdf("../graphcolourings.tvbounds.pdf", width = 10, height = 5)
print(g_tvbounds)
dev.off()
