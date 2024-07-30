#' Global distribution of large submarine canyons
#'
#' Global distribution of large submarine canyons from Harris and Whiteway (2011)
#'
#' @details
#' This dataset has compiled 5849 separate large submarine canyons in the world ocean.
#' @docType data
#' @keywords datasets
#' @format A SpatialLinesDataFrame
#' @source \url{https://doi.org/10.1016/j.margeo.2011.05.008}
#' @references Harris, P.T., Whiteway, T., 2011. Global distribution of large submarine canyons:
#' Geomorphic differences between active and passive continental margins. Marine Geology 285, 69–86.
#' @name canyon
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(sf)
#'
#' bathy <- etopo2022 %>% as.data.frame(xy = TRUE) %>% na.omit
#' ggplot(bathy) +
#'   geom_raster(aes(x=x, y=y, fill=-layer))+
#'   geom_polygon(data=arg, aes(x=X, y=Y, group=PID), fill="bisque2", colour="transparent")+
#'   geom_sf(data=as(eez, "sf"), fill="transparent", colour="red")+
#'   geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=-200, linetype=2, colour="gray50")+
#'   geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=-4000, linetype=1, colour="gray50")+
#'   geom_sf(data=as(canyon, "sf"), size=0.8)+
#'   scale_fill_gradientn(colours=terrain.colors(7))+
#'   scale_x_continuous(expand = expansion(mult = 0))+
#'   scale_y_continuous(expand = expansion(mult = 0))+
#'   labs(x=NULL, y=NULL, fill="Depth\n(m)")+
#'   theme_bw() %+replace% theme(legend.position = "right", legend.key.width =  unit(0.5, 'cm'))
NULL
