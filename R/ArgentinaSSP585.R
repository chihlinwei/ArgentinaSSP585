#' Data package for seafloor climage change (ArgentinaSSP126)
#'
#' ArgentinaSSP126 package contains ensemble model average of CMIP6 historical and ssp126 climate change projections on the seafloor for Argentina EEZ.
#' Yearly mean were calculated from GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-LR, CNRM-ESM2-1, and NorESM2-MM within the Coupled Models Intercomparison Project Phase 6 (cmip6), respectively.
#' The ensemble averages of various model statistics were calculated among the five models for the periods of 1950 to 2000 (historial), 2041 to 2060 (ssp126) and 2081 to 2100 (ssp126).
#' The export POC flux at seafloor was computed from the export production at 100 m (epc100) using the Martin curve (Martin et al., 1987)
#' following the quation: \eqn{Flux = epc100*(depth/export depth)^-0.858}. The depth use \code{\link{etopo2002}} and export depth was
#' set to 100 m. The aragonite (aragsat) and calcite staturation states (calcsat) were calculated as carbonate concentration (co3) divided by
#' carbonate concentration in equilibrium with aragonite (co3satarg) and carbonate concentration in equlibrium with calcite (co3satcalc), respectively.
#' It should be noted that the co3satarg and aragsat were only available from GFDL-ESM4 nd NorESM2-MM.
#' All CMIP6 data were download from \url{https://esgf-node.llnl.gov/search/esgf-llnl/}.
#' ...
#'
#' @author Chih-Lin Wei <chihlinwei@@gmail.com>
#' @docType data
#' @name ArgentinaSSP126
NULL

