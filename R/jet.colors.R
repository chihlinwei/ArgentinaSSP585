#' 'jet.colors' function
#'
#' This function is a wrapper of color interpolation  \code{\link{colorRampPalette}},
#' which takes an integer argument (the required number of colors)
#' and returns a character vector of colors
#'
#' @param n An integer argument (the required number of colors)
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#' @describeIn jet.colors Reversed 'jet.colors'
jet.col.rev <-colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
#' @describeIn jet.colors Easier-on-eye version of 'jet.colors'
jet.colors2 <- colorRampPalette(rev(c("#7F0000", "red2", "red", "orange", "yellow2", "yellow", "#7FFF7F", "cyan", "blue")), bias = 1)
#' @describeIn jet.colors same jet.colors' with "white" color in the middle
jet.colors3 <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","white", "yellow", "#FF7F00", "red", "#7F0000"))
#' @describeIn jet.colors same jet.colors' starts with "white" color
jet.colors4 <-colorRampPalette(c("#FFFFFF", "#7FBFFF", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#' @describeIn jet.colors Same jet.colors in log scale
#' @return a character vector of interpolated colors
#' @author Chih-Lin Wei <chihlinwei@@gmail.com>
#' @examples
#' jet.colors(10)
#' @export
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#' @export
jet.col.rev <-colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))
#' @export
jet.colors2 <- colorRampPalette(rev(c("#7F0000", "red2", "red", "orange", "yellow2", "yellow", "#7FFF7F", "cyan", "blue")), bias = 1)
#' @export
jet.colors3 <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","white", "yellow", "#FF7F00", "red", "#7F0000"))
#' @export
jet.colors4 <-colorRampPalette(c("#FFFFFF", "#7FBFFF", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
