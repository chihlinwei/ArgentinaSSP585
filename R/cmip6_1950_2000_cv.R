#' Coefficient Variation of the Mean Seafloor Climatic Data during 1950 to 2000 among models
#'
#' Coefficient variation of the decadal mean seafloor (bottom-most) projections among models.
#'
#' @details
#' Yearly mean from 1950 to 2000 were calculated from each models and then the decadal mean was calculated from the yearly mean.
#' Coefficient variation was calculated from the decadal mean across all models.
#' \url{https://esgf-node.llnl.gov/search/esgf-llnl/}.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_cv_1950_to_2000}{Export POC flux to seafloor (\%)}
#'   \item{o2_cv_1950_to_2000}{Dissolved oxygen concentration at seafloor (\%)}
#'   \item{ph_cv_1950_to_2000}{pH at seafloor (\%)}
#'   \item{thetao_cv_1950_to_2000}{Potential temperature at seafllor (\%)}
#'   \item{arag_cv_1950_to_2000}{Aragonite Concentration (\%)}
#'   \item{calc_cv_1950_to_2000}{Calcite Concentration (\%)}
#'   \item{co3_cv_1950_to_2000}{Mole Concentration of Carbonate expressed as Carbon in Sea Water (\%)}
#'   \item{co3satarag_cv_1950_to_2000}{Carbonate ion concentration for seawater in equilibrium with pure aragonite (\%)}
#'   \item{co3satcalc_cv_1950_to_2000}{Carbonate ion concentration for seawater in equilibrium with pure calcite (\%)}
#'   \item{aragsat_cv_1950_to_2000}{Aragonite Saturation State (\%)}
#'   \item{calcsat_cv_1950_to_2000}{Calcite Saturation State (\%)}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_1950_2000_cv
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(patchwork)
#' library(tidyr)
#'
#' # Convert raster to data frame and then to list
#' cmip6_list <- subset(cmip6_1950_2000_cv, 1:4) %>%
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
#'     scale_fill_gradientn(colours = jet.colors(10))+
#'     scale_x_continuous(expand = expansion(mult = 0))+
#'     scale_y_continuous(expand = expansion(mult = 0))+
#'     labs(x="Longitude", y="Latitude", fill="")+
#'     coord_fixed(1.52)+
#'     facet_wrap(~ var) +
#'     theme_bw()
#'     })
#'
#' # Wrap ggplot list
#' wrap_plots(gg_list)
NULL
