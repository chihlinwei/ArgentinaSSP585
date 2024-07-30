#' Global distribution of cold-water corals
#'
#' Global distribution of cold-water corals from Freiwald et al. (2017)
#'
#' @details
#' This dataset has compiled occurrence records of cold-water corals in the Argentina EEZ.
#' @docType data
#' @keywords datasets
#' @format A SpatialPointsDataFrame
#' @source \url{http://data.unep-wcmc.org/datasets/3}
#' @references
#' \itemize{
#'   \item{Freiwald A, Rogers A, Hall-Spencer J, Guinotte JM, Davies AJ, Yesson C, Martin CS, Weatherdon LV (2021). Global distribution of cold-water corals (version 5.1).}
#'   \item{Fifth update to the dataset in Freiwald et al. (2004) by UNEP-WCMC, in collaboration with Andre Freiwald and John Guinotte. Cambridge (UK): UN Environment Programme World Conservation Monitoring Centre. Data DOI: \url{https://doi.org/10.34892/72×9-rt61}}
#'   \item{Freiwald A, Fosså JH, Grehan A, Koslow T, Roberts JM (2004). Cold-water coral reefs: out of sight – no longer out of mind. Biodiversity Series 22. Cambridge (UK): UNEP World Conservation Monitoring Centre. 86 pp. URL: \url{https://archive.org/details/coldwatercoralre04frei}}
#'   \item{OSPAR Commission. (2015) OSPAR Threatened and/or Declining Habitats 2015. URL: \url{http://www.ospar.org/work-areas/bdc/species-habitats/list-of-threatened-declining-species-habitats}. Data URL: \url{http://www.emodnet-seabedhabitats.eu/download}.}
#' }
#' @name coral
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
#'   geom_sf(data=as(coral, "sf"), size=0.8)+
#'   scale_fill_gradientn(colours=terrain.colors(7))+
#'   scale_x_continuous(expand = expansion(mult = 0))+
#'   scale_y_continuous(expand = expansion(mult = 0))+
#'   labs(x=NULL, y=NULL, fill="Depth\n(m)")+
#'   theme_bw() %+replace% theme(legend.position = "right", legend.key.width =  unit(0.5, 'cm'))
NULL
