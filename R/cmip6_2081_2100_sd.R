#' Model-Averaging Decadal Standard deviation of Seafloor Climatic Data during 2081 to 2100
#'
#' Emsemble average of the decadal Standard deviation of seafloor (bottom-most) projections from 2081 to 2100.
#' @details
#' Yearly mean projections from 2081 to 2100 were calculated from each models.
#' Decadal standard deviation was calculated from the yearly mean and then averaged across all models.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_sd_2081_to_2100}{Export POC flux to seafloor (mg C m^-2 d^-1)}
#'   \item{o2_sd_2081_to_2100}{Dissolved oxygen concentration at seafloor (mol m^-3)}
#'   \item{ph_sd_2081_to_2100}{pH at seafloor}
#'   \item{thetao_sd_2081_to_2100}{Potential temperature at seafllor (C)}
#'   \item{arag_sd_2081_to_2100}{Aragonite Concentration (mol m-3)}
#'   \item{calc_sd_2081_to_2100}{Calcite Concentration (mol m-3)}
#'   \item{co3_sd_2081_to_2100}{Mole Concentration of Carbonate expressed as Carbon in Sea Water (mol m-3)}
#'   \item{co3satarag_sd_2081_to_2100}{Carbonate ion concentration for seawater in equilibrium with pure aragonite (mol m-3)}
#'   \item{co3satcalc_sd_2081_to_2100}{Carbonate ion concentration for seawater in equilibrium with pure calcite (mol m-3)}
#'   \item{aragsat_sd_2081_to_2100}{Aragonite Saturation State (omega)}
#'   \item{calcsat_sd_2081_to_2100}{Calcite Saturation State (omega)}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_2081_2100_sd
NULL
