#' Model-Averaging of Gradient-based climate velocity magnitude during 2081 to 2100
#'
#' Ensemble average of the seafloor gradient-based climate velocity magnitudes during 2081 to 2100.
#'
#' @details
#' Gradient-based climate velocity magnitude from 2081 to 2100 were calculated from each models and then averaged across all models.
#' Climate velocity was calculated using package gVoCC \url{https://rdrr.io/github/JorGarMol/VoCC/man/gVoCC.html}.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_voccAng_2081_to_2100}{Export POC flux to seafloor (km/yr)}
#'   \item{o2_voccAng_2081_to_2100}{Dissolved oxygen concentration at seafloor (km/yr)}
#'   \item{ph_voccAng_2081_to_2100}{pH at seafloor (km/yr)}
#'   \item{thetao_voccAng_2081_to_2100}{Potential temperature at seafllor (km/yr)}
#'   \item{arag_voccAng_2081_to_2100}{Aragonite Concentration (km/yr)}
#'   \item{calc_voccAng_2081_to_2100}{Calcite Concentration (km/yr)}
#'   \item{co3_voccAng_2081_to_2100}{Mole Concentration of Carbonate expressed as Carbon in Sea Water (km/yr)}
#'   \item{co3satarag_voccAng_2081_to_2100}{Carbonate ion concentration for seawater in equilibrium with pure aragonite (km/yr)}
#'   \item{co3satcalc_voccAng_2081_to_2100}{Carbonate ion concentration for seawater in equilibrium with pure calcite (km/yr)}
#'   \item{aragsat_voccAng_2081_to_2100}{Aragonite saturation State (km/yr)}
#'   \item{calcsat_voccAng_2081_to_2100}{Calcite saturation State (km/yr)}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_2081_2100_voccMeg
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(patchwork)
#' library(tidyr)
#'
#' # Convert raster to data frame and then to list
#' cmip6_list <- subset(cmip6_2081_2100_voccMeg, 1:4) %>%
#'   as.data.frame(xy = TRUE) %>%
#'   gather(-x, -y, key = "var", value = "value") %>%
#'   group_split(var)
#'
#' # Depth
#' bathy <- etopo2022%>% as.data.frame(xy = TRUE) %>% na.omit
#'
#' # ggolot list
#' gg_list = lapply(cmip6_list, function(dat) {
#'   lim <- max(abs(quantile(dat$value, c(0.01, 0.99), na.rm=TRUE)))
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
#'     })
#'
#' # Wrap ggplot list
#' wrap_plots(gg_list)
NULL
