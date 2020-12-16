// Copyright (C) 2018 Dr. Nikolai Knapp, UFZ
//
// This file is part of the MeanShiftR R package.
//
// The MeanShiftR R package is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// MeanShiftR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MeanShiftR If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include <cmath>
#include "LittleFunctionsCollection.h"
using namespace Rcpp;


//' Mean shift clustering
//'
//' @title Mean shift clustering
//' @description
//' Adaptive mean shift clustering to delineate tree crowns from lidar point clouds
//' @param pc Point cloud has to be in matrix format with 3-columns representing X, Y and Z and each row representing one point
//' @param H2CW_fac Factor for the ratio of height to crown width. Determines kernel diameter based on its height above ground.
//' @param H2CL_fac Factor for the ratio of height to crown length. Determines kernel height based on its height above ground.
//' @param UniformKernel Boolean to enable the application of a simple uniform kernel without distance weighting (Default False)
//' @param MaxIter Maximum number of iterations, i.e. steps that the kernel can move for each point. If centroid is not found after all iteration, the last position is assigned as centroid and the processing jumps to the next point
//' @return data.frame with X, Y and Z coordinates of each point in the point cloud and  X, Y and Z coordinates of the centroid to which the point belongs
//' @export
// [[Rcpp::export]]
DataFrame MeanShift_Classical(NumericMatrix pc, double H2CW_fac, double H2CL_fac, bool UniformKernel=false, int MaxIter=20){

  // Create three vectors that all have the length of the incoming point cloud.
  // In these vectors the coodinates of the centroids will be stored
  int nrows = pc.nrow();
  NumericVector centroidx(nrows);
  NumericVector centroidy(nrows);
  NumericVector centroidz(nrows);

  // Loop through all points to process one after the other
  for(int i=0; i<nrows; i++){

    // Initialize variables to store the mean coordinates of all neigbors with the
    // actual coordinates of the focal point from where the kernel starts moving
    double meanx = (double) pc(i,0);
    double meany = (double) pc(i,1);
    double meanz = (double) pc(i,2);

    // Initialize variables to store the old coodinates with unrealistic values of -100
    double oldx = -100.0;
    double oldy = -100.0;
    double oldz = -100.0;

    // Initialize further variables with zero
    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;
    double sump = 0.0;

    double verticalweight = 0.0;
    double horizontalweight = 0.0;
    double weight = 0.0;

    // Keep iterating as long as the centroid (or the maximum number of iterations) is not reached
    int IterCounter = 0;
    do {

      sumx = 0.0;
      sumy = 0.0;
      sumz = 0.0;
      sump = 0.0;

      // Increase the iteration counter
      IterCounter = IterCounter + 1;

      // Calculate cylinder dimensions based on point height
      double r = H2CW_fac * meanz * 0.5;
      double d = H2CW_fac * meanz;
      double h = H2CL_fac * meanz;

      // Remember the coordinate means (kernel position) of previos iteration
       oldx = meanx;
       oldy = meany;
       oldz = meanz;

      // Loop through all points to identify the neighbors of the focal point
      for(int j=0; j<nrows; j++){

        double jx = (double) pc(j,0);
        double jy = (double) pc(j,1);
        double jz = (double) pc(j,2);

        if(InCylinder(jx, jy, jz, r, h, meanx, meany, meanz) == true){

          // If the option of a uniform kernel is set to false calculate the centroid
          // by multiplying all coodinates by their weights, depending on their relative
          // position within the cylinder, summing up the products and dividing by the
          // sum of all weights
          if(UniformKernel == false){
            verticalweight = EpanechnikovFunction(h, meanz, jz);
            horizontalweight = GaussFunction(d, meanx, meany, jx, jy);
            weight = verticalweight * horizontalweight;
            sumx = sumx + weight * jx;
            sumy = sumy + weight * jy;
            sumz = sumz + weight * jz;
            sump = sump + weight;
          }
          // If the option of a uniform kernel is set to true calculate the centroid
          // by summing up all coodinates and dividing by the number of points
          else
          {
            sumx = sumx + jx;
            sumy = sumy + jy;
            sumz = sumz + jz;
            sump = sump + 1.0;
          }
        }
      }
      meanx = sumx / sump;
      meany = sumy / sump;
      meanz = sumz / sump;

    // If the new position equals the previous position (kernel stopped moving), or if the
    // maximum number of iterations is reached, stop the iterations
    }while(meanx != oldx && meany != oldy && meanz != oldz && IterCounter < MaxIter);

    // Store the found position as the centroid position for the focal point
    centroidx[i] = meanx;
    centroidy[i] = meany;
    centroidz[i] = meanz;
  }

 // Return the result as a data.frame with XYZ-coordinates of all points and their corresponding centroids
 return DataFrame::create(_["X"]= pc(_,0),_["Y"]= pc(_,1),_["Z"]= pc(_,2),_["CtrX"]= centroidx,_["CtrY"]= centroidy,_["CtrZ"]= centroidz);
}




















