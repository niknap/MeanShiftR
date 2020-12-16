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

filter_OverlappingCrowns2 <- function(crowns, overlap.threshold=0.9){

  # Package requirements
  # require(rgeos)
  # require(data.table)
  # require(reshape2)

  # Create a boolean vector to label which crowns should be kept
  n.crowns <- nrow(crowns)
  keep.vec <- !logical(n.crowns)
  # Sort crowns by ID
  crowns <- crowns[order(crowns$ID),]
  # Store which pairs of polygons intersect
  intersection.mx <- rgeos::gIntersects(crowns, byid=T)
  # Update the row and column names
  rownames(intersection.mx) <- 1:nrow(crowns)
  colnames(intersection.mx) <- 1:nrow(crowns)
  # Melt the matrix into a table and keep only the intersecting pairs
  intersection.dt <- data.table::data.table(reshape2::melt(intersection.mx))
  setnames(intersection.dt, old=names(intersection.dt), new=c("Crown.i", "Crown.j", "Overlapping"))
  intersection.dt <- subset(intersection.dt, Crown.i != Crown.j & Overlapping == T)
  # Loop over all existing pairwise crown overlaps
  for(my.row in 1:nrow(intersection.dt)){
    #my.row <- 20
    i <- intersection.dt[my.row, Crown.i]
    j <- intersection.dt[my.row, Crown.j]
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
  crowns <- crowns[keep.vec, ]
  return(crowns)
}








