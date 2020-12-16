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



#' Calculate spatial indices for points in 2D space
#'
#' Function that returns spatial indices (plot numbers) for given coordinate pairs.
#' E.g. for hectares res = 100 m
#' @param xcor X-coordinate of the point
#' @param ycor Y-coordinate of the point
#' @param res side length of one spatial subunit
#' @param minx minimum X-coordinate
#' @param miny minimum Y-coordinate
#' @param maxx maximum X-coordinate
#' @param maxy maximum Y-coordinate
#' @param exclude.maxx boolean to set whether the given maxx should be excluded or included in the indexing
#' @return spatial index for the input point
#' @keywords spatial index plot number
#' @author Nikolai Knapp, nikolai.knapp@ufz.de

calc_SpatialIndex <- function(xcor, ycor, res=1, minx=NA, miny=NA, maxx=NA, maxy=NA, exclude.maxx=T){
  # If no coordinate limits are given use the extrema of the data as limits
  if(is.na(minx)){
    minx <- floor(min(xcor, na.rm=T))
  }
  if(is.na(miny)){
    miny <- floor(min(ycor, na.rm=T))
  }
  if(is.na(maxx)){
    maxx <- max(xcor, na.rm=T)
    # Convert absolute to relative maxx
    rel.maxx <- maxx - minx
  } else {
    # Convert absolute to relative maxx
    rel.maxx <- maxx - minx
    if(exclude.maxx==T){
      # Reduce maxx by a tiny number to exclude the upper limit
      rel.maxx <- rel.maxx - 1e-10
    }
  }
  if(is.na(maxy)){
    maxy <- max(ycor, na.rm=T)
  }
  # Convert absolute to relative coordinates
  rel.xcor <- xcor - minx
  rel.ycor <- ycor - miny
  # Count number of cells in X-direction
  # %/%: integer division operator
  nx <- rel.maxx %/% res + 1
  # For each coordinate calculate the cell number in X- and Y-direction
  myx <- rel.xcor %/% res + 1
  myy <- rel.ycor %/% res + 1
  # Calculate the total index
  index <- myx + nx * (myy - 1)
  # Set index values outside the coordinate limits to NA
  index[xcor < minx] <- NA
  index[xcor > maxx] <- NA
  index[ycor < miny] <- NA
  index[ycor > maxy] <- NA
  # Return the result
  return(index)
}







