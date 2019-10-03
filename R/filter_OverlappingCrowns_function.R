#' Remove crown area polygons based on spatial overlap with other polygons
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
#' @param crowns SpatialPolygonsDataFrame of crown projection areas, derived from make_CrownPolygons function
#' @param overlap.threshold Value between 0 and 1 specifying the relative overlap of a pair of crown projection areas (intersection area of two crowns divided by area of a single crown). If the relative overlap is for both crowns larger than this threshold, the crown polygon with the lower tree height value is removed.
#' @return SpatialPolygonsDataFrame with each feature representing the crown projection area of one tree and columns containing various geometric attributes
#' @keywords tree crown projection area polygon overlap cluster artefact error double filter clean remove drop discard
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

filter_OverlappingCrowns <- function(crowns, overlap.threshold){
  require(rgeos)
  n.crowns <- nrow(crowns)
  keep.vec <- !logical(n.crowns)
  # Store which pairs of polygons intersect
  intersection.mx <- gIntersects(crowns, byid=T)
  for(i in 1:n.crowns){
    for(j in 1:n.crowns){
      if(i != j & intersection.mx[i, j] == T){
        # Derive an intersection polygon and calculate its area as a proportion of either of the two crowns
        intersection.polygon <- gIntersection(crowns[i,], crowns[j,], byid=T)
        rel.intersection.i <- gArea(intersection.polygon) / gArea(crowns[i,])
        rel.intersection.j <- gArea(intersection.polygon) / gArea(crowns[j,])
        # Check if threshold is exceeded for both crowns
        if(rel.intersection.i > overlap.threshold & rel.intersection.j > overlap.threshold){
          # Label the lower crown for removal
          head(crowns)
          if(crowns@data[i, "TreeH"] >= crowns@data[j, "TreeH"]){
            keep.vec[j] <- F
          }else{
            keep.vec[i] <- F
          }
        }
      }
    }
  }
  crowns <- crowns[keep.vec, ]
  return(crowns)
}









