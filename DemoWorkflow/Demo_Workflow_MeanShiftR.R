############################################################################################
# Demonstration workflow for individual tree crown delineation from a point cloud using the
# adaptive mean shift 3D algorithm (AMS3D) implemented in MeanShiftR package
############################################################################################

#############
# Preparation
#############

# Load required packages
.libPaths()
require(MeanShiftR)

# Read and inspect example data
require(repmis)
source_data("https://github.com/niknap/MeanShiftR/blob/master/DemoWorkflow/Demo_Data_1ha_Traunstein.rda?raw=True")
head(lid.dt)
head(inv.dt)
nrow(lid.dt)
nrow(inv.dt)

# Display the point cloud (requires function from slidaRtools package available here: https://github.com/niknap/slidaRtools;
# alternatively use the rgl::plot3d function to write your custom 3D plotting code)
require(slidaRtools)
slidaRtools::display.point.cloud.dt(lid.dt, size=2)

####################################
# Data splitting for parallelization
####################################

# Split the point cloud into 25-m tiles with 10-m buffers around each tile
lid.list <- split_BufferedPointCloud(pc.dt=lid.dt, plot.width=25, buffer.width=10)
length(lid.list)

# Check the result by plotting one example tile color coded as core and buffer area
require(slidaRtools)
slidaRtools::display.point.cloud.dt(lid.list[[6]], col.lim=c(0, 1), col.var="Buffer", size=3)

################
# Parallel AMS3D
################

# Apply the parallel AMS3D in the slower "classic" version
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="classic",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))

# Apply the parallel AMS3D in the fast "voxel" version
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="voxel",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))

##################
# Unparallel AMS3D
##################

# Apply the AMS3D in the fast "voxel" version on one large 100-m tile (no parallelization)
lid.list <- split_BufferedPointCloud(pc.dt=lid.dt, plot.width=100, buffer.width=10)
length(lid.list)
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="voxel",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))

# Caution: The execution of the following code may take some time!!!
# Apply the AMS3D in the slow "classic" version on one large 100-m tile (no parallelization)
lid.list <- split_BufferedPointCloud(pc.dt=lid.dt, plot.width=100, buffer.width=10)
length(lid.list)
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="classic",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))

#######################
# Check output visually
#######################

# Count how many tree crown clusters have been found
(Nclus <- length(unique(clus.dt$ID)))

# Display the clustered point cloud with each tree crown in a random color
slidaRtools::display.point.cloud.dt(clus.dt, col.var="ID", col.palette=rainbow(Nclus)[sample(Nclus)], col.lim=c(0, Nclus), size=3)

###############################
# Derive crown projection areas
###############################

# Fit polygons of different shapes for the crown projection areas
ellipses.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="ellipse")
circles.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="circle")
chulls.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="convexhull")

# Create and display a canopy height model (CHM) together with the crown polygons
# (alternative function for CHM generation can be found in the raster package)
require(slidaRtools)
chm.ras <- slidaRtools::raster.from.point.cloud(lid.dt)
plot(chm.ras)
plot(ellipses.spdf, add=T, border="blue")
plot(circles.spdf, add=T, border="red")
plot(chulls.spdf, add=T, border="black")

######################################
# Match crowns with stems in inventory
######################################

# Convert inventory to SpatialPointsDataFrame
inv.spdf <- SpatialPointsDataFrame(coords=cbind(X=inv.dt$X, Y=inv.dt$Y), data=inv.dt)

# Plot inventory data
plot(inv.spdf, cex=inv.spdf@data$DBH, pch=16, col="brown", add=T)

# Add Species column which is required for the following function
inv.spdf@data$Species <- NA

# Find a stem for each crown polygon
ellipses.spdf <- match_CrownsStems(crowns.spdf=ellipses.spdf, stems.spdf=inv.spdf, DBH.min=0.05)
head(ellipses.spdf)

# Plot tree height vs. DBH
plot(ellipses.spdf@data$TreeH ~ ellipses.spdf@data$DBH, xlab="DBH (m)", ylab="Tree height (m)")

# Plot mean tree crown diameter vs. DBH
ellipses.spdf@data$MeanCrownDiameter <- 2 * ellipses.spdf@data$MeanRadius
plot(ellipses.spdf@data$MeanCrownDiameter ~ ellipses.spdf@data$DBH, xlab="DBH (m)", ylab="Crown diameter (m)")






