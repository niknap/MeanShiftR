#' Remove low outlier returns, which originate most probably from tree stems not crowns
#'
#' The AMS3D algorithm can lead to clustering artefacts such that a single tree crown may be
#' assigned to more than one cluster ID and the derived enclosing crown polygons may show very
#' high spatial overlap. This function enables cleaning of the crown polygon dataset by removing
#' redundant polygons which most likely represent just one single tree. Of several
#' overlapping polygons, the polygon with the largest tree height value is kept, while the other
#' polygons are removed.
#' The algorithm involves pairwise comparisons of all polygons which are handed over to the
#' function. Hence, it is recommended to apply it on small tiles (e.g., 1 ha), as processing
#' time may become very slow for large areas.
#' @param pc.dt Point cloud in data.table format containing columns X, Y, Z and ID
#' @param H.quantile Lower end height quantile on which filter factor is applied
#' @param filter.factor which multiplied with the H.quantile height defines the threshold below which a return will be classified as stem.
#' @param rm.class Switch for automatically removing returns belonging to classes. Should be set to "Stem" to remove stem returns, "Crown" to remove crown returns or NA to keep all returns.
#' @return SpatialPolygonsDataFrame with each feature representing the crown projection area of one tree and columns containing various geometric attributes
#' @keywords tree crown projection area polygon overlap cluster artefact error double filter clean remove drop discard
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

filter_StemReturns <- function(pc.dt, H.quantile=0.1, filter.factor=0.9, rm.class=NA){

  # Package requirements
  # require(data.table)

  pc.dt[, LowerQuantileZ := quantile(Z, probs=H.quantile), by=ID]
  pc.dt[, FilterZ := filter.factor * LowerQuantileZ]
  pc.dt[Z >= FilterZ, ReturnClass := "Crown"]
  pc.dt[Z < FilterZ, ReturnClass := "Stem"]
  if(rm.class == "Stem"){
    pc.dt <- subset(pc.dt, ReturnClass != "Stem")
  }else if(rm.class == "Crown"){
    pc.dt <- subset(pc.dt, ReturnClass != "Crown")
  }
  return(pc.dt)
}









