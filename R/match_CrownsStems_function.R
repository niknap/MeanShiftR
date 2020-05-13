#' Match crown projection polygons with tree stem positions
#'
#' The function links remote sensing derived crown projection polygons with inventory derived tree stems.
#' The algorithm sorts all polygons by decreasing area size. Then, starting with the largest crown it checks all
#' stem positions that fall inside the crown projection area and selects the stem with the largest diameter
#' and assigns it to the crown. The stem is removed from the set of available stems and the same procedure
#' is conducted for the next smaller crown and further on until the smallest crown.
#' @param crowns.spdf SpatialPolygonsDataFrame of tree crown projection areas
#' @param stems.spdf SpatialPointsDataFrame of tree stem foot positions (required columns: X, Y, DBH, Species, TreeID)
#' @param DBH.min Minimum stem diameter to be considered
#' @return SpatialPolygonsDataFrame of tree crown projection areas with additional columns containing the attributes of the assigned stem for each crown
#' @keywords tree crown projection area polygon point stem foot position DBH diameter link match inventory ground truth
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

match_CrownsStems <- function(crowns.spdf, stems.spdf, DBH.min=0.05){

  # Package requirements
  # require(data.table)
  # require(rgeos)
  # require(sp)
  # require(raster)

  # Calculate area of each polygon
  Area.vec<- rgeos::gArea(crowns.spdf, byid=T)
  # Sort the crowns for decreasing area
  crowns.spdf <- crowns.spdf[order(Area.vec, decreasing = T),]

  # Filter the inventory data for trees with at least min. DBH
  stems.spdf <- subset(stems.spdf, DBH >= DBH.min)

  # Prepare vectors to store results
  TreeID.vec <- numeric(nrow(crowns.spdf))
  DBH.vec <- numeric(nrow(crowns.spdf))
  Species.vec <- character(nrow(crowns.spdf))
  X.vec <- numeric(nrow(crowns.spdf))
  Y.vec <- numeric(nrow(crowns.spdf))

  # Loop through all crown clusters from large to small and search for the
  # largest stem underneath each crown
  for(i in 1:nrow(crowns.spdf)){
    # i=41
    my.crown.spdf <- crowns.spdf[i, ]
    my.crown.extent <- raster::extent(my.crown.spdf)

    # Find all stems that fall into the crown projection area
    my.pot.stems.spdf <- raster::crop(x=stems.spdf, y=my.crown.extent)
    if(!is.null(my.pot.stems.spdf)){
      over.spdf <- sp::over(x=my.pot.stems.spdf, y=my.crown.spdf)
      my.pot.stems.spdf <- my.pot.stems.spdf[!is.na(over.spdf[, 1]),]
    }

    # Find the largest stem and the associated tree ID
    max.DBH <- NA
    my.ID <- NA
    my.Species <- NA
    my.X <- NA
    my.Y <- NA
    if(!is.null(my.pot.stems.spdf)){
      if(nrow(my.pot.stems.spdf) > 0){
        my.pot.stems.dt <- data.table::data.table(my.pot.stems.spdf@data)
        max.DBH <- max(my.pot.stems.dt$DBH)
        my.ID <- my.pot.stems.dt[DBH == max.DBH, TreeID]
        my.Species <- my.pot.stems.dt[DBH == max.DBH, Species]
        my.X <- my.pot.stems.dt[DBH == max.DBH, X]
        my.Y <- my.pot.stems.dt[DBH == max.DBH, Y]
        # Remove the tree with this ID from the stem table to
        # not use it for multiple crown clusters
        stems.spdf <- subset(stems.spdf, TreeID != my.ID)
      }
    }
    # Store the findings
    TreeID.vec[i] <- my.ID
    DBH.vec[i] <- max.DBH
    Species.vec[i] <- my.Species
    X.vec[i] <- my.X
    Y.vec[i] <- my.Y
  }
  crowns.spdf$StemID <- TreeID.vec
  crowns.spdf$StemX <- X.vec
  crowns.spdf$StemY <- Y.vec
  crowns.spdf$DBH <- DBH.vec
  crowns.spdf$Species <- Species.vec
  return(crowns.spdf)
}













