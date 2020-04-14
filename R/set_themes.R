#'@export
set_theme_chapter1 <- function(){
    library(ggplot2)
    theme_set(theme_void())
    library(gridExtra)
    colors <- c(rgb(.75, 0.75, 0.75), rgb(0.1, 0.5, .8))
    return(list(colors = colors))
}

