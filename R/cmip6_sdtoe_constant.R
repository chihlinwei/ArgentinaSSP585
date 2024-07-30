#' Model-average Constant Time of Emergence based on variability
#'
#' Constant time of emergence was estimated as the first year when future variability constantly
#' (at least 90 percent of the remaining time) exceeding historical variability.
#'
#' @details
#' Constant time of emergence was calculated as the first year when the yearly cumulative future standard deviation
#' constantly (at least 90 percent of the remaining time) exceed the historical yearly standard deviation.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_sdtoe_constant}{Export POC flux to seafloor}
#'   \item{o2_sdtoe_constant}{Dissolved oxygen concentration at seafloor}
#'   \item{ph_sdtoe_constant}{pH at seafloor}
#'   \item{thetao_sdtoe_constant}{Potential temperature at seafllor}
#'   \item{arag_sdtoe_constant}{Aragonite Concentration}
#'   \item{calc_sdtoe_constant}{Calcite Concentration}
#'   \item{co3_sdtoe_constant}{Mole Concentration of Carbonate expressed as Carbon in Sea Water}
#'   \item{co3satarag_sdtoe_constant}{Carbonate ion concentration for seawater in equilibrium with pure aragonite}
#'   \item{co3satcalc_sdtoe_constant}{Carbonate ion concentration for seawater in equilibrium with pure calcite}
#'   \item{aragsat_sdtoe_constant}{Aragonite Saturation State}
#'   \item{calcsat_sdtoe_constant}{Calcite Saturation State}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_sdtoe_constant
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(patchwork)
#' library(tidyr)
#' library(RColorBrewer)
#'
#' # Convert raster to data frame and then to list
#' cmip6_list <- subset(cmip6_sdtoe_constant, 1:4) %>%
#'   as.data.frame(xy = TRUE) %>%
#'   gather(-x, -y, key = "var", value = "value") %>%
#'   group_split(var)
#'
#' # Depth
#' bathy <- etopo2022%>% as.data.frame(xy = TRUE) %>% na.omit
#'
#' # ggolot list
#' gg_list = lapply(cmip6_list, function(dat) {
#'   ggplot(dat) +
#'     geom_raster(aes(x=x, y=y, fill=value))+
#'     geom_polygon(data=arg, aes(x=X, y=Y, group=PID), fill="bisque2", colour="black")+
#'     geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=-200, linetype=2, colour="gray50")+
#'     geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=seq(-6000, -1000, by=1000), linetype=1, colour="gray50")+
#'     scale_fill_gradientn(colours = brewer.pal(10, 'RdYlBu'))+
#'     scale_x_continuous(expand = expansion(mult = 0))+
#'     scale_y_continuous(expand = expansion(mult = 0))+
#'     labs(x="Longitude", y="Latitude", fill="Year")+
#'     coord_fixed(1.52)+
#'     facet_wrap(~ var) +
#'     theme_bw()
#' })
#'
#' # Wrap ggplot list
#' wrap_plots(gg_list)
#'
#' # Overall time of emergence of climate change hazards for export POC flux, dissolved oxygen, pH, and temperature
#' dat <- overlay(subset(cmip6_sdtoe_constant, 1:4), fun=max) %>%
#'   as.data.frame(xy = TRUE)
#'
#' ggplot(dat) +
#'   geom_raster(aes(x=x, y=y, fill=layer))+
#'   geom_polygon(data=arg, aes(x=X, y=Y, group=PID), fill="bisque2", colour="black")+
#'   geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=-200, linetype=2, colour="gray50")+
#'   geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=seq(-6000, -1000, by=1000), linetype=1, colour="gray50")+
#'   scale_fill_gradientn(colours = brewer.pal(10, 'RdYlBu'))+
#'   scale_x_continuous(expand = expansion(mult = 0))+
#'   scale_y_continuous(expand = expansion(mult = 0))+
#'   labs(x="Longitude", y="Latitude", fill="Year")+
#'   coord_fixed(1.52)+
#'   theme_bw()
NULL
