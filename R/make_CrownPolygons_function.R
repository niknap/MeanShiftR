#' Fit polygon hulls around tree crowns in a clustered point cloud to derive crown projection areas
#'
#' The function creates crown projection area polygons from a clustered point cloud of a forest stand.
#' The input point cloud needs to have a column containing tree ID / cluster ID of each point.
#' There are different options for the shape of the polygons.
#' @param pc.dt Point cloud in data.table format containing columns X, Y, Z and ID
#' @param type Shape of the desired polygons ("convexhull", "ellipse" or "circle")
#' @param N.min Minimum number of points in a cluster to be considered as a tree crown
#' @param Ext.min Minimum horizontal extent (in meters, in X- and Y-direction) of a cluster to be considered as a tree crown
#' @param Ext.max Maximum horizontal extent (in meters, in X- and Y-direction) of a cluster to be considered as a tree crown
#' @param proj4string Projection string of class CRS-class
#' @return SpatialPolygonsDataFrame with each feature representing the crown projection area of one tree and columns containing various geometric attributes
#' @keywords tree crown projection area polygons point cloud cluster convex hull ellipse circle perimeter
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

make_CrownPolygons <- function(pc.dt, type="convexhull", N.min=20, Ext.min=2, Ext.max=50, proj4string=CRS(as.character(NA))){

  # type="convexhull"
  # N.min=20
  # Ext.min=2
  # Ext.max=50
  # pc.dt=clus.dt

  # Package requirements
  # require(data.table)
  # require(sp)
  # require(rgeos)
  # require(cluster)
  # require(tripack)
  # require(raster)

  pc.dt <- data.table::data.table(pc.dt)

  # Count the returns per cluster
  pc.dt[, N := .N, by=ID]

  # Subset only clusters that consist of a minimum number of points
  sub.dt <- subset(pc.dt, N >= N.min)

  # Subset only clusters that have a minimum and don't exceed a X- and Y-extent
  sub.dt[, ExtX := max(X, na.rm=T) - min(X, na.rm=T), by=ID]
  sub.dt[, ExtY := max(Y, na.rm=T) - min(Y, na.rm=T), by=ID]
  sub2.dt <- subset(sub.dt, ExtX >= Ext.min & ExtY >= Ext.min & ExtX <= Ext.max & ExtY <= Ext.max)

  # Split into a list of point clouds based on crown IDs
  points.list <- split(sub2.dt, f=sub2.dt$ID)
  # Remove empty elements from the list

  # Pick the selected method
  if(type=="convexhull"){
    # Create a list to store the convex hulls
    hull.list <- list()
    # Loop though the list
    for(i in 1:length(points.list)){
      my.points.dt <- points.list[[i]]
      # Make a SpatialPointsDataFrame from all points
      my.points.spdf <- sp::SpatialPointsDataFrame(coords=cbind(X=my.points.dt$X, Y=my.points.dt$Y), data=my.points.dt, proj4string=proj4string)
      # Collect attributes
      my.ID <- my.points.spdf$ID[1]
      my.CentroidX <- mean(my.points.spdf$X, na.rm=T)
      my.CentroidY <- mean(my.points.spdf$Y, na.rm=T)
      my.CentroidZ <- mean(my.points.spdf$Z, na.rm=T)
      my.TreeH <- max(my.points.spdf$Z, na.rm=T)
      my.CrownBaseH <- min(my.points.spdf$Z, na.rm=T)
      my.CrownLength <- my.TreeH - my.CrownBaseH
      my.NPoints <- my.points.spdf$N[1]
      hull.list[[i]] <- rgeos::gConvexHull(my.points.spdf)
      hull.list[[i]] <- sp::SpatialPolygonsDataFrame(hull.list[[i]], data=data.frame(ID=my.ID, CentroidX=my.CentroidX, CentroidY=my.CentroidY, CentroidZ=my.CentroidZ, TreeH=my.TreeH, CrownBaseH=my.CrownBaseH, CrownLength=my.CrownLength, NPoints=my.NPoints))
    }
    # Bind all list elements together in one SpatialPolygonsDataFrame
    hull.spdf <- do.call(raster::bind, hull.list)
    # Calculate area and perimeter of each polygon
    hull.spdf$ConvexHullArea <- rgeos::gArea(hull.spdf, byid=T)
    hull.spdf$ConvexHullPerimeter <- rgeos::gLength(hull.spdf, byid=T)
    return(hull.spdf)
  } else if(type=="ellipse"){
    # Create a list to store the ellipse hulls
    ellipse.list <- list()
    ps.list <- list()
    # Loop though the list
    for(i in 1:length(points.list)){
      my.points.dt <- points.list[[i]]
      # Make a SpatialPointsDataFrame from all points
      my.points.spdf <- sp::SpatialPointsDataFrame(coords=cbind(X=my.points.dt$X, Y=my.points.dt$Y), data=my.points.dt, proj4string=proj4string)
      # Collect attributes
      my.ID <- my.points.spdf$ID[1]
      my.CentroidX <- mean(my.points.spdf$X, na.rm=T)
      my.CentroidY <- mean(my.points.spdf$Y, na.rm=T)
      my.CentroidZ <- mean(my.points.spdf$Z, na.rm=T)
      my.TreeH <- max(my.points.spdf$Z, na.rm=T)
      my.TreeH <- max(my.points.spdf$Z, na.rm=T)
      my.CrownBaseH <- min(my.points.spdf$Z, na.rm=T)
      my.CrownLength <- my.TreeH - my.CrownBaseH
      my.NPoints <- my.points.spdf$N[1]
      # Combine the coordinates in a matrix
      my.points.mx <- as.matrix(cbind(my.points.spdf$X, my.points.spdf$Y))
      # Fit an ellipse around the data
      ellipse <- cluster::ellipsoidhull(my.points.mx)
      # Calculate border points of the ellipse
      border.points <- cluster::predict.ellipsoid(ellipse)
      # Calculate the center of the ellipse
      center.point <- colMeans(border.points)
      names(center.point) <- c()
      # For each border point calculate the distance to the center
      dist2center <- (rowSums((t(t(border.points)-center.point))^2))^0.5
      # Calculate major and minor semi-axis by finding
      # maximal and minimal distance to the center
      major.semi.axis <- max(dist2center)
      minor.semi.axis <- min(dist2center)
      # Get coordinates of the ellipse point closest and farthest to the center
      min.dist.point <- border.points[dist2center == minor.semi.axis,]
      max.dist.point <- border.points[dist2center == major.semi.axis,]
      # Calculate a mean radius
      mean.radius <- mean(c(major.semi.axis, minor.semi.axis))
      # Collect ellipse attributes
      ellipse.list[[i]] <- data.frame(ID=my.ID, CentroidX=my.CentroidX, CentroidY=my.CentroidY, CentroidZ=my.CentroidZ, TreeH=my.TreeH,
                                      CrownBaseH=my.CrownBaseH, CrownLength=my.CrownLength, NPoints=my.NPoints,
                                      EllipseCenterX=center.point[1], EllipseCenterY=center.point[2], MeanRadius=mean.radius,
                                      MinorSemiAxis=minor.semi.axis, MajorSemiAxis=major.semi.axis)
      # Make a polygon
      p = Polygon(border.points)
      # Add it to the Polygons list
      ps.list[[i]] <- Polygons(list(p), ID=i)
    }
    # Convert polygons to spatial polygons
    sps <- sp::SpatialPolygons(ps.list)
    # Bind all list elements together in one data.table
    ellipse.dt <- do.call(raster::bind, ellipse.list)
    # Convert them to a spatial polygon dataframe
    ellipse.spdf <- sp::SpatialPolygonsDataFrame(sps, data=ellipse.dt)
    # Calculate area and perimeter of each polygon
    ellipse.spdf$EllipseArea <- rgeos::gArea(ellipse.spdf, byid=T)
    ellipse.spdf$EllipsePerimeter <- rgeos::gLength(ellipse.spdf, byid=T)
    return(ellipse.spdf)
  } else if(type=="circle"){
    # Create a list to store the circle hulls
    circle.list <- list()
    ps.list <- list()
    # Loop though the list
    for(i in 1:length(points.list)){
      my.points.dt <- points.list[[i]]
      # Collect attributes
      my.ID <- points.list[[i]]$ID[1]
      my.CentroidX <- mean(my.points.dt$X, na.rm=T)
      my.CentroidY <- mean(my.points.dt$Y, na.rm=T)
      my.CentroidZ <- mean(my.points.dt$Z, na.rm=T)
      my.TreeH <- max(my.points.dt$Z, na.rm=T)
      my.CrownBaseH <- min(my.points.dt$Z, na.rm=T)
      my.CrownLength <- my.TreeH - my.CrownBaseH
      my.NPoints <- my.points.dt$N[1]
      # Add some very small random noise to the coordinates to avoid errors with duplicate
      # coordinate combinations in the circumcircle function
      my.points.dt[, X := X + 0.01*(runif(nrow(my.points.dt))-0.5)]
      my.points.dt[, Y := Y + 0.01*(runif(nrow(my.points.dt))-0.5)]
      # Find the smallest circumcircle of the points
      circle <- tripack::circumcircle(x=my.points.dt$X, y=my.points.dt$Y)
      # Extract center and radius of the circle
      CircleCenterX <- circle$x
      CircleCenterY <- circle$y
      CircleRadius <- circle$radius
      # Make the center a spatial object
      CircleCenter.sp <- sp::SpatialPoints(coords=data.frame(X=CircleCenterX, Y=CircleCenterY), proj4string=proj4string)
      # Create the circle polygon
      ps.list[[i]] <- rgeos::gBuffer(CircleCenter.sp, width=CircleRadius, quadsegs=15)
      circle.list[[i]] <- data.frame(ID=my.ID, CentroidX=my.CentroidX, CentroidY=my.CentroidY, CentroidZ=my.CentroidZ, TreeH=my.TreeH,
                                     CrownBaseH=my.CrownBaseH, CrownLength=my.CrownLength, NPoints=my.NPoints,
                                     CircleCenterX=CircleCenterX, CircleCenterY=CircleCenterY, Radius=CircleRadius)
    }
    # Bind all list elements together in one data.table
    circle.dt <- do.call(raster::bind, circle.list)
    # Bind all list elements together in one SpatialPolygonsDataFrame.
    circle.spdf <- do.call(raster::bind, ps.list)
    # Convert them to a spatial polygon dataframe
    circle.spdf <- sp::SpatialPolygonsDataFrame(circle.spdf, data=circle.dt)
    # Calculate area and perimeter of each polygon
    circle.spdf$CircleArea <- rgeos::gArea(circle.spdf, byid=T)
    circle.spdf$CirclePerimeter <- rgeos::gLength(circle.spdf, byid=T)
    return(circle.spdf)
  }
}



















