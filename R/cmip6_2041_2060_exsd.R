#' Exposure to Climate Change Hazard by 2041 to 2060
#'
#' Magnitude of climate change by 2041 to 2060 in the unit of historical variability
#'
#' @details
#' Exposure to climate change hazard was calculated as \code{\link{cmip6_2041_2060_ch}} divided by \code{\link{cmip6_1950_2000_sd}}.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_exsd_2041_to_2060}{Exposure to climate change hazard in export POC flux to seafloor}
#'   \item{o2_exsd_2041_to_2060}{Exposure to climate change hazard in dissolved oxygen concentration at seafloor}
#'   \item{ph_exsd_2041_to_2060}{Exposure to climate change hazard in pH at seafloor}
#'   \item{thetao_exsd_2041_to_2060}{Exposure to climate change hazard in potential temperature at seafllor}
#'   \item{arag_exsd_2041_to_2060}{Exposure to climate change hazard in aragonite Concentration}
#'   \item{calc_exsd_2041_to_2060}{Exposure to climate change hazard in calcite Concentration}
#'   \item{co3_exsd_2041_to_2060}{Exposure to climate change hazard in mole Concentration of Carbonate expressed as Carbon in Sea Water}
#'   \item{co3satarag_exsd_2041_to_2060}{Exposure to climate change hazard in carbonate ion concentration for seawater in equilibrium with pure aragonite}
#'   \item{co3satcalc_exsd_2041_to_2060}{Exposure to climate change hazard in carbonate ion concentration for seawater in equilibrium with pure calcite}
#'   \item{aragsat_exsd_2041_to_2060}{Exposure to climate change hazard in aragonite Saturation State}
#'   \item{calcsat_exsd_2041_to_2060}{Exposure to climate change hazard in calcite Saturation State}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_2041_2060_exsd
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(patchwork)
#' library(tidyr)
#'
#' # Exposure to Climate Hazards
#'
#' # Convert raster to data frame and then to list
#' cmip6_list <- subset(cmip6_2041_2060_exsd, 1:4) %>%
#'   as.data.frame(xy = TRUE) %>%
#'   gather(-x, -y, key = "var", value = "value") %>%
#'   group_split(var)
#'
#' # Depth
#' bathy <- etopo2022%>% as.data.frame(xy = TRUE) %>% na.omit
#'
#' # ggolot list
#' gg_list = lapply(cmip6_list, function(dat) {
#'   lim <- max(abs(quantile(dat$value, c(0.0001, 0.9999), na.rm=TRUE)))
#'   ggplot(dat) +
#'     geom_raster(aes(x=x, y=y, fill=value))+
#'     geom_polygon(data=arg, aes(x=X, y=Y, group=PID), fill="bisque2", colour="black")+
#'     geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=-200, linetype=2, colour="gray50")+
#'     geom_contour(data=bathy, aes(x=x, y=y, z=layer), breaks=seq(-6000, -1000, by=1000), linetype=1, colour="gray50")+
#'     scale_fill_gradientn(colours = jet.colors3(10), limits=c(-lim, lim))+
#'     scale_x_continuous(expand = expansion(mult = 0))+
#'     scale_y_continuous(expand = expansion(mult = 0))+
#'     labs(x="Longitude", y="Latitude", fill="")+
#'     coord_fixed(1.52)+
#'     facet_wrap(~ var) +
#'     theme_bw()
#' })
#'
#' # Wrap ggplot list
#' wrap_plots(gg_list)
NULL
