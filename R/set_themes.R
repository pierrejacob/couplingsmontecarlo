#'@export
set_theme_chapter1 <- function(){
    library(ggplot2)
    library(gridExtra)
    theme_set(theme_void())
    colors <- c(rgb(.75, 0.75, 0.75), rgb(0.1, 0.5, .8))
    return(list(colors = colors))
}

#'@export
set_theme_chapter2 <- function(){
    library(ggplot2)
    library(gridExtra)
    theme_set(theme_void())
    colors <- c(rgb(.1, 0.5, 0.3), rgb(0.9, 0.8, 0.1))
    return(list(colors = colors))
}
