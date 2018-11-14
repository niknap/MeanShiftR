#' Split point cloud into smaller point clouds with buffer area around
#'
#' The function splits one large point cloud into several smaller point clouds.
#' It allows to specify a buffer width around the core area where points are included.
#' @param pc.dt Point cloud in data.table format containing columns X, Y and Z
#' @param plot.width Width of the core area in meters
#' @param buffer.width Width of the buffer around the core area in meters
#' @return List of sub point clouds and buffer and core points labeled in boolean column "Buffer"
#' @keywords point cloud split buffer area plot subset parallel
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

split_BufferedPointCloud <- function(pc.dt, plot.width, buffer.width){

  # Package requirements
  require(data.table)
  require(plyr)

  # Convert to data.table
  pc.dt <- data.table(pc.dt)

  # Calculate the absolute lower left corner of the point cloud
  abs.llX <- round_any(min(pc.dt$X, na.rm=T), accuracy=plot.width, f=floor)
  abs.llY <- round_any(min(pc.dt$Y, na.rm=T), accuracy=plot.width, f=floor)

  # Calculate spatial indices for small subplots
  if(is.data.table(pc.dt)){pc.dt[, sBPC_SpatID := calc_SpatialIndex(xcor=X, ycor=Y, res=plot.width, minx=abs.llX, miny=abs.llY)]}

  # Calculate coordinates of the lower left plot corners
  pc.dt[, sBPC_llX := round_any(X, accuracy=plot.width, f=floor)]
  pc.dt[, sBPC_llY := round_any(Y, accuracy=plot.width, f=floor)]

  # Make a unique data.table with all coordinate corner combinations and spatial indices
  coord.dt <- unique(subset(pc.dt, select=c("sBPC_SpatID", "sBPC_llX", "sBPC_llY")))
  setorderv(coord.dt, cols=c("sBPC_SpatID", "sBPC_llX", "sBPC_llY"))
  setnames(coord.dt, old=names(coord.dt), new=c("sBPC_nSpatID", "sBPC_nllX", "sBPC_nllY"))

  # Create a copy of the point cloud which will become the result and to which the buffers will be added successively
  result.dt <- copy(pc.dt)
  result.dt[, Buffer := 0]

  # Add the lower buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX]
  buf.dt[, sBPC_nllY := sBPC_llY+plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, Y < sBPC_nllY & Y >= sBPC_nllY - buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the lower left buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX+plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY+plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer.width & Y < sBPC_nllY & Y >= sBPC_nllY - buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the left buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX+plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper left buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX+plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY-plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X < sBPC_nllX & X >= sBPC_nllX - buffer.width & Y > sBPC_nllY + plot.width & Y <= sBPC_nllY + plot.width + buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX]
  buf.dt[, sBPC_nllY := sBPC_llY-plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, Y > sBPC_nllY + plot.width & Y <= sBPC_nllY + plot.width + buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the upper right buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX-plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY-plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + plot.width & X <= sBPC_nllX + plot.width + buffer.width & Y > sBPC_nllY + plot.width & Y <= sBPC_nllY + plot.width + buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the right buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX-plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + plot.width & X <= sBPC_nllX + plot.width + buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Add the lower right buffer points
  buf.dt <- copy(pc.dt)
  # Calculate the neighbor plot's lower left coordinates
  buf.dt[, sBPC_nllX := sBPC_llX-plot.width]
  buf.dt[, sBPC_nllY := sBPC_llY+plot.width]
  # Add the ID of the neigbor plot
  buf.dt <- merge(buf.dt, coord.dt, by=c("sBPC_nllX", "sBPC_nllY"), all.x=T)
  # Subset points inside the buffer
  buf.dt <- subset(buf.dt, X > sBPC_nllX + plot.width & X <= sBPC_nllX + plot.width + buffer.width & Y < sBPC_nllY & Y >= sBPC_nllY - buffer.width)
  # Add the buffer points to the result
  buf.dt[, Buffer := 1]
  buf.dt[, sBPC_SpatID := sBPC_nSpatID]
  buf.dt <- subset(buf.dt, select=names(result.dt))
  result.dt <- rbind(result.dt, buf.dt)
  rm(buf.dt)

  # Remove all rows with NA as spatial index
  result.dt <- subset(result.dt, !is.na(sBPC_SpatID))

  # Split into a list of separate point clouds with buffers based on spatial index
  setorderv(result.dt, cols=c("sBPC_SpatID", "Buffer", "X", "Y"))
  result.list <- split(result.dt, by="sBPC_SpatID")

  return(result.list)

}




