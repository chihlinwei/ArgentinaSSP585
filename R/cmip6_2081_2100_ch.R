#' Seafloor Climate Change by 2081 to 2100
#'
#' Emsemble average of the decadal mean seafloor (bottom-most) projections during 2081 to 2100
#' minus the decadal mean projections during 1950 to 2000.
#'
#' @details
#' Climate change (difference of decadal mean) between the period of 1950 to 2000 and the period 2041 to 2080 were calculated
#' for each models. Emsemble average of climate changes was calculated across all models.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_av_2081_to_2100_minus_1950_to_2000}{Export POC flux to seafloor (mg C m^-2 d^-1)}
#'   \item{o2_av_2081_to_2100_minus_1950_to_2000}{Dissolved oxygen concentration at seafloor (mol m^-3)}
#'   \item{ph_av_2081_to_2100_minus_1950_to_2000}{Difference in pH at seafloor}
#'   \item{thetao_av_2081_to_2100_minus_1950_to_2000}{Potential temperature at seafllor (C)}
#'   \item{arag_av_2081_to_2100_minus_1950_to_2000}{Aragonite concentration (mol m-3)}
#'   \item{calc_av_2081_to_2100_minus_1950_to_2000}{Calcite concentration (mol m-3)}
#'   \item{co3_av_2081_to_2100_minus_1950_to_2000}{Mole concentration of Carbonate expressed as Carbon in Sea Water (mol m-3)}
#'   \item{co3satarag_av_2081_to_2100_minus_1950_to_2000}{Carbonate ion concentration for seawater in equilibrium with pure aragonite (mol m-3)}
#'   \item{co3satcalc_av_2081_to_2100_minus_1950_to_2000}{Carbonate ion concentration for seawater in equilibrium with pure calcite (mol m-3)}
#'   \item{aragsat_av_2081_to_2100_minus_1950_to_2000}{Aragonite Saturation State (omega)}
#'   \item{calcsat_av_2081_to_2100_minus_1950_to_2000}{Calcite Saturation State (omega)}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_2081_2100_ch
NULL
