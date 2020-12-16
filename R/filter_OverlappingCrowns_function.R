# Copyright (C) 2018 Dr. Nikolai Knapp, UFZ
#
# This file is part of the MeanShiftR R package.
#
# The MeanShiftR R package is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# MeanShiftR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MeanShiftR If not, see <http://www.gnu.org/licenses/>.



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

  # Package requirements
  # require(rgeos)

  # Create a boolean vector to label which crowns should be kept
  n.crowns <- nrow(crowns)
  keep.vec <- !logical(n.crowns)
  # Store which pairs of polygons intersect
  intersection.mx <- rgeos::gIntersects(crowns, byid=T)
  for(i in 1:n.crowns){
    for(j in 1:n.crowns){
      if(i != j & intersection.mx[i, j] == T){
        # Derive an intersection polygon and calculate its area as a proportion of either of the two crowns
        intersection.polygon <- rgeos::gIntersection(crowns[i,], crowns[j,], byid=T)
        rel.intersection.i <- rgeos::gArea(intersection.polygon) / rgeos::gArea(crowns[i,])
        rel.intersection.j <- rgeos::gArea(intersection.polygon) / rgeos::gArea(crowns[j,])
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









