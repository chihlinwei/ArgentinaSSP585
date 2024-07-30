#' Global distribution of large seamounts
#'
#' Global distribution of large seamounts from Kim and Wessel (2011)
#'
#' @details
#' This dataset has compiled 24643  separate large seamounts in the world ocean.
#' @docType data
#' @keywords datasets
#' @format A SpatialPointsDataFrame
#' @source \url{https://doi.org/10.1111/j.1365-246X.2011.05076.x}
#' @references Kim, S.-S., Wessel, P., 2011. New global seamount census from altimetry-derived gravity data.
#' Geophysical Journal International 186, 615–631.
#' @name seamount
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
#'   geom_sf(data=as(seamount, "sf"), size=0.8)+
#'   scale_fill_gradientn(colours=terrain.colors(7))+
#'   scale_x_continuous(expand = expansion(mult = 0))+
#'   scale_y_continuous(expand = expansion(mult = 0))+
#'   labs(x=NULL, y=NULL, fill="Depth\n(m)")+
#'   theme_bw() %+replace% theme(legend.position = "right", legend.key.width =  unit(0.5, 'cm'))
NULL
