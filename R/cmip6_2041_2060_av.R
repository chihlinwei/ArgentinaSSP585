#' Model-Averaging of Decadal Mean Seafloor Climatic Data during 2041 to 2060
#'
#' Ensemble average of the decadal mean seafloor (bottom-most) projections from 2041 to 2060.
#'
#' @details
#' Yearly mean projections from 2041 to 2060 were calculated from each models
#' Decadal mean was calculated from the yearly mean and than averaged across all models.
#' @docType data
#' @keywords datasets
#' @format A RasterBrick object of 11 raster layers:
#' \describe{
#'   \item{epc_av_2041_to_2060}{Export POC flux to seafloor (mg C m^-2 d^-1)}
#'   \item{o2_av_2041_to_2060}{Dissolved oxygen concentration at seafloor (mol m^-3)}
#'   \item{ph_av_2041_to_2060}{pH at seafloor}
#'   \item{thetao_av_2041_to_2060}{Potential temperature at seafllor (C)}
#'   \item{arag_av_2041_to_2060}{Aragonite Concentration (mol m-3)}
#'   \item{calc_av_2041_to_2060}{Calcite Concentration (mol m-3)}
#'   \item{co3_av_2041_to_2060}{Mole Concentration of Carbonate expressed as Carbon in Sea Water (mol m-3)}
#'   \item{co3satarag_av_2041_to_2060}{Carbonate ion concentration for seawater in equilibrium with pure aragonite (mol m-3)}
#'   \item{co3satcalc_av_2041_to_2060}{Carbonate ion concentration for seawater in equilibrium with pure calcite (mol m-3)}
#'   \item{aragsat_av_2041_to_2060}{Aragonite Saturation State (omega)}
#'   \item{calcsat_av_2041_to_2060}{Calcite Saturation State (omega)}
#' }
#' @source \url{https://esgf-node.llnl.gov/search/esgf-llnl/}
#' @name cmip6_2041_2060_av
NULL
