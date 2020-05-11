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
# load igraph package to deal with graphs
library(igraph)
# function to plot graph
custom_plot_graph <- function(g) plot(g, layout=layout_with_kk, vertex.label = NA)

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
all_colours <- colorRampPalette(RColorBrewer::brewer.pal(q, name = "Paired"))(q)
# all_colours[1:2] <- graphsettings$colors
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


## Markov chain
single_kernel <- function(g){
    # choose a vertex
    ivertex <- sample(1:nvertices, 1)
    n_i <- neighbors(g, ivertex)
    # get colors not in neighbors
    legal_colours <- setdiff(all_colours, V(g)$color[n_i])
    V(g)$color[ivertex] <- sample(x = legal_colours, size = 1)
    return(list(g = g, ivertex = ivertex))
}

## initialize chains
g <- rinit(g)
ghistory <- list()
vertexhistory <- list()
## advance one chain for 'lag' steps
nmcmc <- 5e2
for (imcmc in 1:nmcmc){
    res_ <- single_kernel(g)
    g <- res_$g
    ghistory[[imcmc]] <- g
    vertexhistory[[imcmc]] <- res_$ivertex
}


# pdf("../traceplot.graphcolourings.pdf")
par(mfrow = c(4,4), mar = c(1,1,1,1))
for (i in 1:16){
    plot(ghistory[[i]], layout=layout_with_kk, vertex.label = NA, mark.groups = vertexhistory[[i]],
         mark.col = rgb(0,0,0,0.5), mark.border = rgb(0,0,0,1))
    title(paste0("iteration ", i))
}
# dev.off()

