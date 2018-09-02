############################################################################################
# Demonstration workflow for individual tree crown delineation from a point cloud using the
# adaptive mean shift 3D algorithm (AMS3D) implemented in MeanShiftR package
############################################################################################

# Load required packages
.libPaths()
require(MeanShiftR)

# Set working directory
# wd <- "C:\\Users\\knappn\\owncloud\\Laptop\\RScripts\\PackageDevelopment_MeanShiftR\\TestFiles\\"
# setwd(wd)

# Read and inspect example data
lid.dt <- readRDS("Traunstein.rectangle.lidar.dt.rds")
display.point.cloud.dt(lid.dt)
nrow(lid.dt)/180000

# Subset the data for testing
lid.dt <- subset(lid.dt, X < 400 & Y < 200)
display.point.cloud.dt(lid.dt, size=2)

# Split the point cloud into 50-m tiles with 10-m buffers around
lid.list <- split_BufferedPointCloud(pc.dt=lid.dt, plot.width=100, buffer.width=10)
length(lid.list)

# Check the result by plotting one example tile color coded as core and buffer area
require(slidaRtools)
display.point.cloud.dt(lid.list[[10]], col.lim=c(0, 1), col.var="Buffer", size=3)

# Apply the parallel AMS3D in the fast "voxel" version
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="voxel",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))
260/18
55/18

# Apply the parallel AMS3D in the slower "classic" version
system.time(clus.dt <- parallel_MeanShift(pc.list=lid.list, lib.path=.libPaths()[1], frac.cores=0.5, version="classic",
                                          H2CW=0.3, H2CL=0.4, max.iter=20, buffer.width=10, minz=2, ctr.ac=2))
1200/18

# Count how many tree crown clusters have been found
(Nclus <- length(unique(clus.dt$ID)))

# Display the clustered point cloud with each tree crown in a random color
display.point.cloud.dt(clus.dt, col.var="ID", col.palette=rainbow(Nclus)[sample(Nclus)], col.lim=c(0, Nclus), size=3)

# Fit polygons of different shapes for the crown projection areas
circles.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="circle")
ellipses.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="ellipse")
chulls.spdf <- make_CrownPolygons(pc.dt=clus.dt, type="convexhull")

# Create and display a canopy height model together with the crown polygons
chm.ras <- raster.from.point.cloud(lid.dt)
plot(chm.ras)
plot(circles.spdf, add=T, border="red")
plot(ellipses.spdf, add=T, border="blue")
plot(chulls.spdf, add=T, border="black")

# Read example inventory data
inv.dt <- data.table(readRDS("Traunstein.rectangle.inventory.df.rds"))

# Subset the data for testing
inv.dt <- subset(inv.dt, X < 100 & Y < 100)

# Convert to SpatialPointsDataFrame
inv.spdf <- SpatialPointsDataFrame(coords=cbind(X=inv.dt$X, Y=inv.dt$Y), data=inv.dt)

# Plot inventory data
plot(inv.spdf, cex=inv.spdf$DBH, pch=16, col="brown", add=T)

# Find a stem for each crown polygon
ellipses.spdf <- match_CrownsStems(crowns.spdf=ellipses.spdf, stems.spdf=inv.spdf, DBH.min=0.05)
head(ellipses.spdf)

# Plot tree height vs. DBH
plot(ellipses.spdf$TreeH ~ ellipses.spdf$DBH)

# Plot mean tree crown diameter vs. DBH
ellipses.spdf$MeanCrownDiameter <- 2 * ellipses.spdf$MeanRadius
plot(ellipses.spdf$MeanCrownDiameter ~ ellipses.spdf$DBH)






